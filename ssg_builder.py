import os
import json
import math
import itertools
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
        .pagination { margin-top: 20px; display: flex; gap: 10px; flex-wrap: wrap; }
        .page-link { padding: 5px 10px; border: 1px solid #ccc; text-decoration: none; color: black; margin-bottom: 5px;}
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
                <th>Gene / Feature ID</th>
                <th>Full Locus Signature</th>
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
        <input type="text" id="geneSearch" placeholder="Search for a gene (e.g., TP53) or feature...">
        <ul id="searchResults" class="results-list"></ul>
    </div>

    <h2>Browse Biological Interactions</h2>
    {% for super_cat, sub_cats in directory_tree.items() %}
        <h3>📁 {{ super_cat }}</h3>
        <ul>
        {% for sub_cat, chroms in sub_cats.items() %}
            <li>📁 <strong>{{ sub_cat }}</strong> 
                ({% for chrom in chroms %}<a href="{{ super_cat|replace(' ', '_')|replace('/', '-') }}/{{ sub_cat|replace(' ', '_')|replace('/', '-') }}/{{ chrom }}_page_1.html">{{ chrom }}</a>{% if not loop.last %}, {% endif %}{% endfor %})
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
                        li.innerHTML = `<a href="${url}">View occurrences of ${gene}</a>`;
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

def get_main_super_cat(f1, f2):
    coding_kws = ['cds', 'exon', 'gene', 'transcript', 'protein', 'start_codon', 'stop_codon', 'utr', 'selenocysteine']
    c1_coding = any(k in f1.lower() for k in coding_kws)
    c2_coding = any(k in f2.lower() for k in coding_kws)
    
    if c1_coding and c2_coding: return "Protein-Coding / Protein-Coding"
    elif c1_coding or c2_coding: return "Protein-Coding / Non-Coding"
    else: return "Non-Coding / Non-Coding"

def route_overlap(sig):
    items = [x.strip() for x in str(sig).split('|')]
    
    # Filter out structural fluff if specific features exist
    generic = {'gene', 'transcript', 'exon'}
    specific = [x for x in items if x.lower() not in generic]
    
    if not specific:
        if 'exon' in [x.lower() for x in items]: specific = ['exon']
        elif 'transcript' in [x.lower() for x in items]: specific = ['transcript']
        else: specific = ['gene']
        
    seen = set()
    cleaned = [x.capitalize() for x in specific if not (x.lower() in seen or seen.add(x.lower()))]
    
    routes = []
    
    if len(cleaned) > 2:
        # Route 1: The full complex string to the Complex Interactions folder
        complex_sub = " vs ".join(cleaned)[:60]
        routes.append(("Complex Interactions", complex_sub))
        
        # Route 2: Extract all clean pairs and map them to the 3 main folders
        for combo in itertools.combinations(cleaned, 2):
            super_cat = get_main_super_cat(combo[0], combo[1])
            sub_cat = f"{combo[0]} vs {combo[1]}"
            routes.append((super_cat, sub_cat))
            
    elif len(cleaned) == 2:
        super_cat = get_main_super_cat(cleaned[0], cleaned[1])
        sub_cat = f"{cleaned[0]} vs {cleaned[1]}"
        routes.append((super_cat, sub_cat))
        
    else: # Length is exactly 1 (e.g. CDS vs CDS)
        f = cleaned[0]
        super_cat = get_main_super_cat(f, f)
        sub_cat = f"{f} vs {f}"
        routes.append((super_cat, sub_cat))
        
    return routes

def get_f1(row):
    return row['gene_ids'] if row['gene_ids'] != 'none' else row['feature_ids']

def generate_static_site_chunked(csv_path, output_dir="genomic-overlap-atlas"):
    base_path = Path(output_dir)
    base_path.mkdir(exist_ok=True)
    
    search_index = defaultdict(set)
    directory_tree = defaultdict(lambda: defaultdict(list))
    
    print(f"📦 Pass 1: Chunking and Multi-Mapping dataset...")
    chunk_size = 250000
    for chunk_idx, chunk in enumerate(pd.read_csv(csv_path, chunksize=chunk_size, sep=None, engine='python')):
        print(f"   -> Processing rows {chunk_idx * chunk_size} to {(chunk_idx + 1) * chunk_size}...")
        
        # Apply the multi-mapping array to each row
        chunk['routes'] = chunk['feature_signature'].apply(route_overlap)
        chunk['f1'] = chunk.apply(get_f1, axis=1)
        
        # Explode duplicates the row for every route it matches (UX Magic)
        exploded_chunk = chunk.explode('routes')
        exploded_chunk['super_cat'] = exploded_chunk['routes'].apply(lambda x: x[0])
        exploded_chunk['sub_cat'] = exploded_chunk['routes'].apply(lambda x: x[1])
        
        for (super_c, sub_c, chrom), group in exploded_chunk.groupby(['super_cat', 'sub_cat', 'chromosome']):
            safe_super = super_c.replace(" ", "_").replace("/", "-")
            safe_sub = sub_c.replace(" ", "_").replace("/", "-")
            
            dir_path = base_path / safe_super / safe_sub
            dir_path.mkdir(parents=True, exist_ok=True)
            
            out_df = pd.DataFrame({
                'chr': group['chromosome'],
                'start': group['start'],
                'end': group['end'],
                'feature1': group['f1'],
                'feature2': group['feature_signature'], # Keep the full string so they can see the complexity
                'length': group['length']
            })
            
            tsv_file = dir_path / f"{chrom}_raw_data.tsv"
            out_df.to_csv(tsv_file, mode='a', header=not tsv_file.exists(), sep='\t', index=False)
            
            if chrom not in directory_tree[super_c][sub_c]:
                directory_tree[super_c][sub_c].append(chrom)
                
            target_url = f"{safe_super}/{safe_sub}/{chrom}_page_1.html"
            genes = group[group['gene_ids'] != 'none']['gene_ids'].unique()
            for g_str in genes:
                for g in str(g_str).split(','):
                    search_index[g.strip().upper()].add(target_url)

    print("🗂️ Pass 2: Sorting lengths and generating HTML pages...")
    table_tpl = Template(TABLE_TEMPLATE_HTML)
    index_tpl = Template(INDEX_TEMPLATE_HTML)
    ROWS_PER_PAGE = 5000
    total_processed = 0

    ordered_tree = {
        "Protein-Coding / Protein-Coding": directory_tree.get("Protein-Coding / Protein-Coding", {}),
        "Protein-Coding / Non-Coding": directory_tree.get("Protein-Coding / Non-Coding", {}),
        "Non-Coding / Non-Coding": directory_tree.get("Non-Coding / Non-Coding", {}),
        "Complex Interactions": directory_tree.get("Complex Interactions", {})
    }
    
    ordered_tree = {k: v for k, v in ordered_tree.items() if v}

    for super_c, sub_cats in ordered_tree.items():
        safe_super = super_c.replace(" ", "_").replace("/", "-")
        for sub_c, chroms in sub_cats.items():
            safe_sub = sub_c.replace(" ", "_").replace("/", "-")
            
            chroms.sort(key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else 999)
            
            for chrom in chroms:
                dir_path = base_path / safe_super / safe_sub
                tsv_file = dir_path / f"{chrom}_raw_data.tsv"
                
                df_chrom = pd.read_csv(tsv_file, sep='\t')
                df_chrom = df_chrom.sort_values(by='length', ascending=False)
                df_chrom.to_csv(tsv_file, sep='\t', index=False)
                
                rows = df_chrom.to_dict('records')
                total_rows = len(rows)
                total_processed += total_rows
                total_pages = math.ceil(total_rows / ROWS_PER_PAGE)
                
                for page_num in range(1, total_pages + 1):
                    start_idx = (page_num - 1) * ROWS_PER_PAGE
                    page_rows = rows[start_idx : start_idx + ROWS_PER_PAGE]
                    
                    html = table_tpl.render(
                        super_cat=super_c, sub_cat=sub_c, chrom=chrom,
                        rows=page_rows, current_page=page_num, total_pages=total_pages,
                        total_rows=total_rows, tsv_filename=tsv_file.name
                    )
                    with open(dir_path / f"{chrom}_page_{page_num}.html", "w", encoding="utf-8") as f:
                        f.write(html)

    with open(base_path / "search_index.json", "w") as f:
        json.dump({k: list(v) for k, v in search_index.items()}, f)

    with open(base_path / "index.html", "w", encoding="utf-8") as f:
        f.write(index_tpl.render(directory_tree=ordered_tree))
        
    print(f"✅ Clean SSG Build Complete! Multi-mapped {total_processed} overlaps to 4 strict folders.")

if __name__ == "__main__":
    generate_static_site_chunked("genomic_overlap_segments.csv")
