"""Build docs/Ilha_study_questions.pdf: a 10-question self-quiz (answers at the back)
plus the 20 hard questions from 04_hard_questions.md. Run with the project venv:
    PATH=../.venv/bin:$PATH python3 build_study_pdf.py
"""
import html, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CSS = os.path.join(HERE, "..", "ilha_erw_model", "report", "doc.css")
TOPDF = os.path.join(HERE, "..", "ilha_erw_model", "report", "topdf.py")

QUIZ = [
    ("The rock has a CO2 potential of about 320 kg per tonne. Why is that not what we expect from the trial?",
     "It is a ceiling: what the rock could capture if every grain dissolved and every cation left as bicarbonate, which takes decades. "
     "358 t x 320 kg is about 115 t. In 30 months only about 7% dissolves (about 8 t) and only about one eighth of that leaves as bicarbonate (about 1 t gross). "
     "Trucking, spreading, rock carbon and sulfur cost about 2.4 t, so the trial is about -1.3 t net."),
    ("What share of the rock's cations dissolves by month 30 in the central case, and what is the honest range?",
     "About 7% (central). Range about 3% (pessimistic) to 13% (optimistic). 7% is the middle case, not the optimistic one."),
    ("At soil pH 5.4, roughly what share of the released cations leaves as bicarbonate? And at pH 6.2?",
     "About 10% at pH 5.4; about 41% at pH 6.2 (22% at 5.8). The rest neutralises soil acidity. That is why rising pH speeds up carbon removal later."),
    ("If the December screening finds soil pH 5.0 instead of 5.4, does the trial's carbon get better or worse? Why?",
     "Worse. The rock dissolves a bit faster (8.1% vs 6.7% at 40 t/ha) but only about 4% leaves as bicarbonate, so the trial net goes from about -1.3 t to about -1.7 t. "
     "More acid means more dissolving but less carbon. The farm benefit (pH rise) is bigger."),
    ("Which lab result tells us most about how much of the rock dissolves in 30 months?",
     "The grain-size test (PSD) of the Tecnovia fines, including laser diffraction of the part finer than 63 micrometres. "
     "Those grains are about 12% of the mass but give about 70% of the early dissolution, like icing sugar versus a sugar cube."),
    ("What is the field factor, why does it matter, and where will we get it from?",
     "How fast rock weathers in a real field compared with a lab beaker. It is the single biggest uncertainty (rank correlation 0.83). "
     "We will estimate it from the mesocosms (soil columns in the lab) and from the first field rounds, fitted back through the model."),
    ("Name three measurements that should show a rate effect within 30 months, and two that will not.",
     "Will: soil pH, exchangeable Ca + Mg, total Ca + Mg (soil cores), lysimeter alkalinity under 40 and 80 t/ha, mesocosm leachate. "
     "Will not: the dissolved fraction from soil mass balance (detectable only above about 0.17 vs 0.07 expected), the Ti and Zr tracers, and the streams W1 to W6 (signal at most 0.01 meq/L)."),
    ("Why 16 plots rather than 12?",
     "At medium variability 16 plots detect the dose response with power about 0.95, 12 plots about 0.8. "
     "Caveat: both figures assume the central field factor; if the rock weathers slower, power drops. Drop to 12 only if the screening shows low variability."),
    ("Why is '40 t/ha versus control' a weak test, and what do we use as the main result instead?",
     "Comparing single rates uses only 8 plots and needs a Holm correction for several comparisons: power about 0.35 even with 16 plots. "
     "The main result is the regression slope of the response on rate across all 16 plots, with blocks."),
    ("Are the trace metals a problem?",
     "Not from the basalt so far, but open. With placeholder soil values the natural background already exceeds the Portuguese sludge-law screening limits for Ni, Cr, Cu and Zn, which is common in volcanic soils. "
     "80 t/ha adds about +10% Ni and +20% Cr. Which limits apply to rock dust on volcanic soil must be settled with Juan and the regional authority after the December screening."),
]


def esc(s):
    return html.escape(s)


def md_inline(s):
    s = esc(s)
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"\*(.+?)\*", r"<i>\1</i>", s)
    return s


def hard_questions():
    txt = open(os.path.join(HERE, "04_hard_questions.md")).read()
    out = []
    for line in txt.splitlines():
        if line.startswith("## "):
            out.append(f"<h3>{esc(line[3:])}</h3>")
        elif re.match(r"^\*\*\d+\.", line):
            out.append(f"<p class='q'>{md_inline(line)}</p>")
        elif line.strip() and not line.startswith("# ") and not line.startswith("Numbers are from"):
            out.append(f"<p>{md_inline(line)}</p>")
    return "\n".join(out)


def build():
    q = "\n".join(f"<li><p class='q'>{esc(a)}</p><div class='lines'></div></li>" for a, _ in QUIZ)
    a = "\n".join(f"<li><p><b>{esc(qq)}</b></p><p>{esc(aa)}</p></li>" for qq, aa in QUIZ)
    page = f"""<!doctype html><html><head><meta charset="utf-8"><title>Ilha study questions</title>
<link rel="stylesheet" href="{CSS}">
<style>
.q {{ font-weight: 700; margin-top: 6pt; }}
.lines {{ height: 30mm; border-bottom: 0.4pt solid #999; background: repeating-linear-gradient(#fff 0 7.4mm, #bbb 7.4mm 7.6mm); margin-bottom: 4pt; }}
ol li {{ page-break-inside: avoid; }}
.pb {{ page-break-before: always; }}
</style></head><body>
<div class="titleblock"><div class="doclabel">Ilha Carbon · Lomba da Maia basalt rate trial · Study questions</div>
<h1>Questions to work through before the team meeting</h1>
<div class="sub">19 November 2026. Part A: a 10-question self-quiz with space to write (answers in Part C). Part B: the 20 hardest questions the team is likely to ask, with short answers.</div></div>
<div class="box thick"><p><b>Placeholder values: replace with screening and lab data.</b> Every number here comes from the model run on literature soil values and the Marques Britas basalt. None is a site result. Say "placeholder" whenever you quote one.</p></div>
<h2>Part A. Self-quiz</h2>
<p>Answer in your own words first, then check Part C.</p>
<ol>{q}</ol>
<h2 class="pb">Part B. The 20 hard questions, with short answers</h2>
{hard_questions()}
<h2 class="pb">Part C. Answers to the self-quiz</h2>
<ol>{a}</ol>
</body></html>"""
    src = os.path.join(HERE, "Ilha_study_questions.html")
    out = os.path.join(HERE, "Ilha_study_questions.pdf")
    open(src, "w").write(page)
    subprocess.run([sys.executable, TOPDF, src, out, "Ilha Carbon · study questions · placeholder values"], check=True)
    print(out)


if __name__ == "__main__":
    build()
