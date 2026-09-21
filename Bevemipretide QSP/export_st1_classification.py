# -*- coding: utf-8 -*-
import docx, json
path = r"C:\Users\ssr92\OneDrive\Desktop\My Projects\Bevemipretide PD QSP\New publication plots\Manuscript_updated.docx"
d = docx.Document(path)
t = d.tables[7]
by_class = {}
for row in t.rows[1:]:
    cells = [c.text.strip() for c in row.cells]
    name, cls = cells[0], cells[5]
    by_class.setdefault(cls, []).append(name)

counts = {k: len(v) for k, v in by_class.items()}
counts["Total"] = sum(counts.values())

out = {"counts": counts, "params_by_class": by_class}
with open("ST1_classification.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

print("Counts:", counts)
print("Saved ST1_classification.json")
