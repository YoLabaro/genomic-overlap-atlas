import os
import math
import pandas as pd

INPUT_FILES = {
    "internal": "strict_internal_feature_overlaps.tsv",
    "dual": "ultimate_dual_function_overlaps.tsv",
    "triads": "enhancer_promoter_ctcf_triads.tsv"
}

OUTPUT_DIR = "genomic_overlap_atlas_site"
ROWS_PER_PAGE = 5000

def get_ucsc_link(chrom, start, end):
    span = end - start
    pad = max(100, int(span * 0.1))
    u_start = max(0, start - pad)
    u_end = end + pad
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}:{u_start}-{u_end}"

def generate_html_table(df_chunk, title, breadcrumbs, current_page, total_pages):
    rows_html = ""
    for _, row in df_chunk.iterrows():
        chrom = str(row.get('chrom', row.get('Chromosome', 'chr1')))
        start = int(row.get('overlap_start', row.get('start', 0)))
        end = int(row.get('overlap_end', row.get('end', start + int(row.get('length', 0)))))
        length = int(row.get('length', 0))
        ucsc_url = get_ucsc_link(chrom, start, end)
        
        feat_id = str(row.get('id_1', row.get('feature_id', row.get('feature_1_id', ''))))
        feat_id_b = str(row.get('id_2', row.get('feature_id_b', row.get('feature_2_id', ''))))
        full_id_str = f"{feat_id} | {feat_id_b}" if feat_id_b and feat_id_b != 'nan' else feat_id
        context = str(row.get('biological_context', row.get('feature_combination', row.get('feature_type', ''))))
        locus_str = f"{chrom}:{start:,}-{end:,}"

        rows_html += f"""
        <tr>
            <td class="mono">{locus_str}</td>
            <td><strong>{full_id_str}</strong></td>
            <td><span class="context-tag">{context}</span></td>
            <td><span class="length-badge">{length:,} bp</span></td>
            <td><a href="{ucsc_url}" target="_blank" class="btn-ucsc">🧬 View in UCSC</a></td>
        </tr>
        """

    pagination_html = ""
    if total_pages > 1:
        pagination_html = '<div class="pagination">'
        for p in range(1, total_pages + 1):
            if p == current_page:
                pagination_html += f'<span class="page-btn active">{p}</span>'
            else:
                p_name = "index.html" if p == 1 else f"page_{p}.html"
                pagination_html += f'<a href="{p_name}" class="page-btn">{p}</a>'
        pagination_html += '</div>'

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{title} - Genomic Overlap Atlas</title>
    <style>
        :root {{
            --bg-color: #f8fafc;
            --card-bg: #ffffff;
            --text-main: #1e293b;
            --text-muted: #64748b;
            --border-color: #e2e8f0;
            --primary: #2563eb;
            --primary-hover: #1d4ed8;
            --accent-bg: #eff6ff;
        }}
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 0; padding: 30px; background: var(--bg-color); color: var(--text-main); }}
        .container {{ max-width: 1400px; margin: auto; background: var(--card-bg); padding: 35px; border-radius: 12px; box-shadow: 0 4px 20px rgba(0,0,0,0.05); border: 1px solid var(--border-color); }}
        .breadcrumbs {{ font-size: 14px; color: var(--text-muted); margin-bottom: 15px; }}
        .breadcrumbs a {{ color: var(--primary); text-decoration: none; font-weight: 500; }}
        .breadcrumbs a:hover {{ text-decoration: underline; }}
        h1 {{ margin-top: 0; font-size: 26px; color: #0f172a; }}
        .toolbar {{ display: flex; justify-content: space-between; align-items: center; margin: 25px 0 15px 0; flex-wrap: wrap; gap: 15px; }}
        .search-box {{ padding: 10px 14px; width: 320px; border: 1px solid var(--border-color); border-radius: 6px; font-size: 14px; outline: none; transition: border-color 0.2s; }}
        .search-box:focus {{ border-color: var(--primary); box-shadow: 0 0 0 3px rgba(37,99,235,0.1); }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
        th, td {{ padding: 14px 16px; border-bottom: 1px solid var(--border-color); text-align: left; font-size: 14px; }}
        th {{ background-color: #f1f5f9; font-weight: 600; color: #334155; position: sticky; top: 0; }}
        tr:hover {{ background-color: #f8fafc; }}
        .mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; color: #475569; }}
        .length-badge {{ background: var(--accent-bg); color: var(--primary); padding: 4px 10px; border-radius: 20px; font-weight: 600; font-size: 13px; display: inline-block; }}
        .context-tag {{ background: #f1f5f9; color: #475569; padding: 4px 8px; border-radius: 4px; font-size: 13px; }}
        .btn-ucsc {{ background: #eff6ff; color: var(--primary); padding: 7px 14px; border-radius: 6px; text-decoration: none; font-weight: 600; font-size: 13px; display: inline-block; border: 1px solid #bfdbfe; transition: all 0.2s; }}
        .btn-ucsc:hover {{ background: var(--primary); color: white; border-color: var(--primary); }}
        .pagination {{ margin-top: 30px; display: flex; gap: 6px; flex-wrap: wrap; align-items: center; }}
        .page-btn {{ padding: 8px 14px; border: 1px solid var(--border-color); border-radius: 6px; text-decoration: none; color: var(--primary); font-size: 14px; background: #fff; font-weight: 500; }}
        .page-btn.active {{ background: var(--primary); color: #fff; border-color: var(--primary); }}
        .page-btn:hover:not(.active) {{ background: #f1f5f9; }}
        .meta-info {{ font-size: 14px; color: var(--text-muted); }}
    </style>
</head>
<body>
    <div class="container">
        <div class="breadcrumbs">{breadcrumbs}</div>
        <h1>{title}</h1>
        <div class="toolbar">
            <div class="meta-info">Sorted by absolute genomic length (longest to shortest). Page {current_page} of {total_pages}.</div>
            <input type="text" id="tableSearch" class="search-box" placeholder="🔍 Live filter table (e.g. gene, locus)..." onkeyup="filterTable()">
        </div>
        <table>
            <thead>
                <tr>
                    <th>Locus</th>
                    <th>Gene / Feature ID</th>
                    <th>Biological Context</th>
                    <th>Overlap Size</th>
                    <th>Action</th>
                </tr>
            </thead>
            <tbody id="dataTable">
                {rows_html}
            </tbody>
        </table>
        {pagination_html}
    </div>

    <script>
    function filterTable() {{
        let input = document.getElementById("tableSearch");
        let filter = input.value.toLowerCase();
        let table = document.getElementById("dataTable");
        let tr = table.getElementsByTagName("tr");
        for (let i = 0; i < tr.length; i++) {{
            let txtValue = tr[i].textContent || tr[i].innerText;
            if (txtValue.toLowerCase().indexOf(filter) > -1) {{
                tr[i].style.display = "";
            }} else {{
                tr[i].style.display = "none";
            }}
        }}
    }}
    </script>
</body>
</html>
"""

def main():
    print("🚀 Building Styled Atlas 2.0...")
    
    subfolder_mapping = {
        "internal": "internal",
        "dual": "dual",
        "triads": "higher_complexity"
    }

    _super_categories = [
        "proteincoding-proteincoding",
        "proteincoding-nonproteincoding",
        "nonproteincoding-nonproteincoding"
    ]

    for key, filename in INPUT_FILES.items():
        if os.path.exists(filename):
            print(f"📖 Processing {filename}...")
            df = pd.read_csv(filename, sep='\t')
            
            if 'Chromosome' in df.columns and 'chrom' not in df.columns:
                df['chrom'] = df['Chromosome']
            if 'overlap_start' not in df.columns and 'start' in df.columns:
                df['overlap_start'] = df['start']
            if 'overlap_end' not in df.columns and 'end' in df.columns:
                df['overlap_end'] = df['end']

            if 'overlap_start' in df.columns and 'overlap_end' in df.columns:
                df['length'] = df['overlap_end'] - df['overlap_start']
            elif 'triple_overlap_bp' in df.columns:
                df['length'] = df['triple_overlap_bp']
                df['overlap_start'] = df.get('overlap_start', 0)
                df['overlap_end'] = df['overlap_start'] + df['length']
            else:
                df['length'] = 0
                
            df = df.sort_values(by='length', ascending=False).reset_index(drop=True)
            sub_dir_name = subfolder_mapping[key]
            total_rows = len(df)
            total_pages = math.ceil(total_rows / ROWS_PER_PAGE)

            for super_cat in _super_categories:
                target_dir = os.path.join(super_cat, sub_dir_name)
                os.makedirs(target_dir, exist_ok=True)
                
                title = f"{super_cat.replace('-', ' / ').title()} — {sub_dir_name.replace('_', ' ').title()}"
                breadcrumbs = f'<a href="../../../index.html">Home</a> &gt; <span>{super_cat}</span> &gt; <span>{sub_dir_name}</span>'
                
                for page_num in range(1, total_pages + 1):
                    start_idx = (page_num - 1) * ROWS_PER_PAGE
                    end_idx = start_idx + ROWS_PER_PAGE
                    chunk = df.iloc[start_idx:end_idx]
                    
                    filename_html = "index.html" if page_num == 1 else f"page_{page_num}.html"
                    out_path = os.path.join(target_dir, filename_html)
                    
                    html_content = generate_html_table(chunk, title, breadcrumbs, page_num, total_pages)
                    with open(out_path, "w") as f:
                        f.write(html_content)

    root_index = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Genomic Overlap Atlas 2.0</title>
    <style>
        :root {
            --bg-color: #f8fafc;
            --card-bg: #ffffff;
            --text-main: #1e293b;
            --text-muted: #64748b;
            --border-color: #e2e8f0;
            --primary: #2563eb;
        }
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 0; padding: 50px 20px; background: var(--bg-color); color: var(--text-main); }
        .container { max-width: 1000px; margin: auto; background: var(--card-bg); padding: 50px; border-radius: 16px; box-shadow: 0 10px 30px rgba(0,0,0,0.06); border: 1px solid var(--border-color); }
        h1 { color: #0f172a; margin-top: 0; font-size: 32px; display: flex; align-items: center; gap: 12px; }
        p.subtitle { color: var(--text-muted); font-size: 16px; line-height: 1.6; margin-bottom: 40px; }
        .category-grid { display: grid; grid-template-columns: 1fr; gap: 24px; }
        .card { background: #ffffff; border: 1px solid var(--border-color); border-radius: 12px; padding: 30px; transition: all 0.25s ease; box-shadow: 0 2px 4px rgba(0,0,0,0.02); }
        .card:hover { transform: translateY(-3px); box-shadow: 0 12px 24px rgba(37,99,235,0.08); border-color: #93c5fd; }
        .card h2 { margin-top: 0; color: #0f172a; font-size: 20px; display: flex; align-items: center; gap: 10px; }
        .card p { color: var(--text-muted); font-size: 14px; line-height: 1.5; margin-bottom: 20px; }
        .links { display: flex; gap: 12px; flex-wrap: wrap; }
        .btn { background: #eff6ff; color: var(--primary); padding: 10px 16px; border-radius: 8px; text-decoration: none; font-weight: 600; font-size: 13px; border: 1px solid #bfdbfe; transition: all 0.2s; }
        .btn:hover { background: var(--primary); color: white; border-color: var(--primary); }
        .btn-complex { background: #fefce8; color: #ca8a04; border-color: #fde047; }
        .btn-complex:hover { background: #ca8a04; color: white; border-color: #ca8a04; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Genomic Overlap Atlas 2.0</h1>
        <p class="subtitle">An advanced interactive exploration portal for high-confidence structural genomic feature interactions, dual-function overlaps, and regulatory triads across the human genome (hg38). Features live filtering and direct gene-anchored UCSC Genome Browser links.</p>
        
        <div class="category-grid">
            <div class="card">
                <h2>📁 Protein-Coding / Protein-Coding</h2>
                <p>Explore gene fusions, alternative splicing conflicts, read-through transcripts, and dense exonic stacking.</p>
                <div class="links">
                    <a href="proteincoding-proteincoding/internal/index.html" class="btn">Internal Overlaps</a>
                    <a href="proteincoding-proteincoding/dual/index.html" class="btn">Dual-Function</a>
                    <a href="proteincoding-proteincoding/higher_complexity/index.html" class="btn btn-complex">Higher Complexity (Triads)</a>
                </div>
            </div>

            <div class="card">
                <h2>📁 Protein-Coding / Non-Coding</h2>
                <p>Examine regulatory control mechanisms, silencers, enhancers, and structural repeats interrupting coding regions.</p>
                <div class="links">
                    <a href="proteincoding-nonproteincoding/internal/index.html" class="btn">Internal Overlaps</a>
                    <a href="proteincoding-nonproteincoding/dual/index.html" class="btn">Dual-Function</a>
                    <a href="proteincoding-nonproteincoding/higher_complexity/index.html" class="btn btn-complex">Higher Complexity (Triads)</a>
                </div>
            </div>

            <div class="card">
                <h2>📁 Non-Coding / Non-Coding</h2>
                <p>Investigate gene deserts, pure regulatory networks, long non-coding RNA interactions, and genome organization.</p>
                <div class="links">
                    <a href="nonproteincoding-nonproteincoding/internal/index.html" class="btn">Internal Overlaps</a>
                    <a href="nonproteincoding-nonproteincoding/dual/index.html" class="btn">Dual-Function</a>
                    <a href="nonproteincoding-nonproteincoding/higher_complexity/index.html" class="btn btn-complex">Higher Complexity (Triads)</a>
                </div>
            </div>
        </div>
    </div>
</body>
</html>
"""
    with open("index.html", "w") as f:
        f.write(root_index)
    print("✨ Atlas 2.0 generation complete!")

if __name__ == "__main__":
    main()
