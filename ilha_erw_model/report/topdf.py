import sys, asyncio
from playwright.sync_api import sync_playwright
src, out, label = sys.argv[1], sys.argv[2], sys.argv[3]
landscape = len(sys.argv) > 4 and sys.argv[4] == 'land'
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page()
    import pathlib
    pg.goto(pathlib.Path(src).resolve().as_uri()); pg.wait_for_timeout(600)
    foot = f'<div style="font-family:Inter,Arial;font-size:7.5px;width:100%;padding:0 18mm;display:flex;justify-content:space-between;color:#444"><span>{label}</span><span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>'
    pg.pdf(path=out, format='A4', landscape=landscape, print_background=True, display_header_footer=True,
           header_template='<div></div>', footer_template=foot, prefer_css_page_size=True)
    b.close()
