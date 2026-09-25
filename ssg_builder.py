import os
import json
import math
import pandas as pd
from pathlib import Path
from collections import defaultdict
from jinja2 import Template

TABLE_TEMPLATE_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{{ super_cat }} | {{ sub_cat }} | {{ chrom }}</title>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .header { display: flex; justify-content: space-between; align-items: center; }
        .btn { padding: 8px 16px; background-color: #0073e6; color: white; text-decoration: none; border-radius: 4px; }
        .pagination { margin-top: 20px; display: flex; gap: 10px; }
        .page-link { padding: 5px 10px; border: 1px solid #ccc; text-decoration: none; color: black; }
        .page-link.active { background-color: #0073e6; color: white; border-color: #0073e6; }
    </style>
</head>
<body>
    <div class="header">
        <h1>{{ sub_cat }} ({{ chrom }})</h1>
        <div>
            <a href="../../index.html" class="btn">🏠 Home</a>
            <a href="{{ tsv_filename }}" class="btn" download>💾 Download Full TSV</a>
        </div>
    </div>
    <p>Showing page {{ current_page }} of {{ total_pages }} ({{ total_rows }} total overlaps found).</p>
    
    <table id="overlapTable" class="display">
        <thead>
            <tr>
                <th>Locus</th>
                <th>Feature 1</th>
                <th>Feature 2</th>
                <th>Overlap Size (bp)</th>
                <th>Action</th>
            </tr>
        </thead>
        <tbody>
            {% for row in rows %}
            <tr>
                <td>{{ row.chr }}:{{ row.start }}-{{ row.end }}</td>
                <td>{{ row.feature1 }}</td>
                <td>{{ row.feature2 }}</td>
                <td>{{ row.length }}</td>
                <td><a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={{ row.chr }}:{{ row.start }}-{{ row.end }}" target="_blank">🧬 View in UCSC</a></td>
            </tr>
            {% endfor %}
        </tbody>
    </table>

    <div class="pagination">
        {% for p in range(1, total_pages + 1) %}
            <a href="page_{{ p }}.html" class="page-link {% if p == current_page %}active{% endif %}">Page {{ p }}</a>
        {% endfor %}
    </div>

    <script>
        $(document).ready(function() {$('#overlapTable').DataTable({
                "order": [[ 3, "desc" ]] 
            });
        });
    </script>
</body>
</html>
"""

INDEX_TEMPLATE_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Genomic Overlap Atlas</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; max-width: 900px; margin: auto; }
        .search-container { margin-bottom: 30px; }
        input[type="text"] { width: 100%; padding: 12px; font-size: 16px; border: 1px solid #ccc; border-radius: 4px; }
        .folder-list { list-style: none; padding: 0; }
        .folder-list li { margin: 10px 0; font-size: 18px; }
        .results-list { margin-top: 10px; background: #f9f9f9; border: 1px solid #eee; padding: 10px; display: none; }
    </style>
</head>
<body>
    <h1>Genomic Overlap Atlas</h1>
    
    <div class="search-container">
        <input type="text" id="geneSearch" placeholder="Search for a gene (e.g., TP53)...">
        <ul id="searchResults" class="results-list"></ul>
    </div>

    <h2>Browse Biological Interactions</h2>
    {% for super_cat, sub_cats in directory_tree.items() %}
        <h3>📁 {{ super_cat }}</h3>
        <ul>
        {% for sub_cat, chroms in sub_cats.items() %}
            <li>📁 <strong>{{ sub_cat }}</strong> 
                ({% for chrom in chroms %}<a href="{{ super_cat|replace(' ', '_') }}/{{ sub_cat|replace(' ', '_') }}/{{ chrom }}_page_1.html">{{ chrom }}</a>{% if not loop.last %}, {% endif %}{% endfor %})
            </li>
        {% endfor %}
        </ul>
    {% endfor %}

    <script>
        let searchIndex = {};
        fetch('search_index.json')
            .then(response => response.json())
            .then(data => searchIndex = data);

        document.getElementById('geneSearch').addEventListener('input', function(e) {
            const query = e.target.value.toUpperCase();
            const resultsBox = document.getElementById('searchResults');
            resultsBox.innerHTML = '';
            
            if (query.length < 2) {
                resultsBox.style.display = 'none';
                return;
            }

            let found = false;
            for (const [gene, urls] of Object.entries(searchIndex)) {
                if (gene.includes(query)) {
                    found = true;
                    urls.forEach(url => {
                        let li = document.createElement('li');
                        li.innerHTML = `<a href="${url}">View overlaps containing ${gene}</a>`;
                        resultsBox.appendChild(li);
                    });
                }
            }
            resultsBox.style.display = found ? 'block' : 'none';
        });
    </script>
</body>
</html>
"""

def generate_static_site_from_csv(csv_path, output_dir="genomic-overlap-atlas"):
    base_path = Path(output_dir)
    base_path.mkdir(exist_ok=True)
    
    print(f"📖 Reading genomic data from {csv_path}...")
    df = pd.read_csv(csv_path, sep=None, engine='python')
    
    grouped_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    search_index = defaultdict(set)
    
    # Flexible column mapping based on standard bedtools/sweep-line headers
    chr_col = next((c for c in df.columns if 'chr' in c.lower() or c == 'chrom'), df.columns[0])
    start_col = next((c for c in df.columns if 'start' in c.lower()), df.columns[1])
    end_col = next((c for c in df.columns if 'end' in c.lower()), df.columns[2])
    f1_col = next((c for c in df.columns if 'name' in c.lower() or 'feature1' in c.lower() or 'gene' in c.lower()), df.columns[3])
    f2_col = next((c for c in df.columns if 'feature2' in c.lower() or 'overlap' in c.lower()), df.columns[4] if len(df.columns) > 4 else f1_col)

    for _, row in df.iterrows():
        chrom = str(row[chr_col])
        try:
            start = int(row[start_col])
            end = int(row[end_col])
        except ValueError:
            continue
            
        f1 = str(row[f1_col])
        f2 = str(row[f2_col])
        length = end - start
        
        # Automatic biological categorization rule
        coding_keywords = ['cds', 'exon', 'gene', 'transcript', 'protein']
        is_f1_coding = any(k in f1.lower() for k in coding_keywords)
        is_f2_coding = any(k in f2.lower() for k in coding_keywords)
        
        if is_f1_coding and is_f2_coding:
            super_cat = "Protein-Coding / Protein-Coding"
        elif is_f1_coding or is_f2_coding:
            super_cat = "Protein-Coding / Non-Coding"
        else:
            super_cat = "Non-Coding / Non-Coding"
            
        sub_cat = f"{f1[:15]} vs {f2[:15]}"
        
        record = {
            'chr': chrom, 'start': start, 'end': end,
            'feature1': f1, 'feature2': f2, 'length': length,
            'super_cat': super_cat, 'sub_cat': sub_cat
        }
        
        grouped_data[super_cat][sub_cat][chrom].append(record)
        
        safe_super = super_cat.replace(" ", "_")
        safe_sub = sub_cat.replace(" ", "_")
        target_url = f"{safe_super}/{safe_sub}/{chrom}_page_1.html"
        search_index[f1.upper()].add(target_url)

    table_tpl = Template(TABLE_TEMPLATE_HTML)
    index_tpl = Template(INDEX_TEMPLATE_HTML)
    
    directory_tree_for_index = defaultdict(lambda: defaultdict(list))
    ROWS_PER_PAGE = 5000
    total_processed = 0

    for super_c, sub_cats in grouped_data.items():
        safe_super = super_c.replace(" ", "_")
        for sub_c, chroms in sub_cats.items():
            safe_sub = sub_c.replace(" ", "_")
            for chrom, rows in chroms.items():
                directory_tree_for_index[super_c][sub_c].append(chrom)
                
                dir_path = base_path / safe_super / safe_sub
                dir_path.mkdir(parents=True, exist_ok=True)
                
                rows.sort(key=lambda x: x['length'], reverse=True)
                total_processed += len(rows)
                
                tsv_filename = f"{chrom}_raw_data.tsv"
                pd.DataFrame(rows).to_csv(dir_path / tsv_filename, sep='\t', index=False)
                
                total_rows = len(rows)
                total_pages = math.ceil(total_rows / ROWS_PER_PAGE)
                
                for page_num in range(1, total_pages + 1):
                    start_idx = (page_num - 1) * ROWS_PER_PAGE
                    end_idx = start_idx + ROWS_PER_PAGE
                    page_rows = rows[start_idx:end_idx]
                    
                    html_content = table_tpl.render(
                        super_cat=super_c, sub_cat=sub_c, chrom=chrom,
                        rows=page_rows, current_page=page_num, total_pages=total_pages,
                        total_rows=total_rows, tsv_filename=tsv_filename
                    )
                    
                    with open(dir_path / f"{chrom}_page_{page_num}.html", "w", encoding="utf-8") as f:
                        f.write(html_content)

    with open(base_path / "search_index.json", "w") as f:
        json.dump({k: list(v) for k, v in search_index.items()}, f)

    with open(base_path / "index.html", "w", encoding="utf-8") as f:
        f.write(index_tpl.render(directory_tree=directory_tree_for_index))
        
    print(f"✅ SSG Build Complete! Successfully processed and hosted {total_processed} real genomic overlaps.")

if __name__ == "__main__":
    generate_static_site_from_csv("genomic_pairwise_overlaps.csv")
