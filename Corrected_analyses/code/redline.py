"""Minimal tracked-change (w:ins / w:del) editing helpers for WordprocessingML paragraphs and tables."""
import copy
import difflib
import re
import zipfile
from lxml import etree as E

WNS = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
NS = {'w': WNS}
XML_SPACE = '{http://www.w3.org/XML/1998/namespace}space'


def q(tag): return '{%s}%s' % (WNS, tag)


class Tracker:
    def __init__(self, author, date, start_id=90000):
        self.author, self.date, self.next_id = author, date, start_id

    def mark(self, tag):
        e = E.Element(q(tag))
        self.next_id += 1
        e.set(q('id'), str(self.next_id)); e.set(q('author'), self.author); e.set(q('date'), self.date)
        return e

    # ------------------------------------------------------------ run plumbing
    @staticmethod
    def normalize(p):
        """Give every run a single content child so text runs can be split cleanly."""
        for r in list(p.iter(q('r'))):
            kids = [c for c in r if c.tag != q('rPr')]
            if len(kids) <= 1: continue
            rpr = r.find(q('rPr')); parent = r.getparent(); pos = parent.index(r)
            parent.remove(r)
            for k, child in enumerate(kids):
                nr = E.Element(q('r'))
                if rpr is not None: nr.append(copy.deepcopy(rpr))
                nr.append(child); parent.insert(pos + k, nr)

    @staticmethod
    def segments(p):
        out, pos = [], 0
        for r in p.iter(q('r')):
            if any(a.tag == q('del') for a in r.iterancestors()): continue
            t = r.find(q('t'))
            if t is None: continue
            n = len(t.text or ''); out.append((r, pos, pos + n)); pos += n
        return out

    @classmethod
    def text(cls, p):
        return ''.join(r.find(q('t')).text or '' for r, _, _ in cls.segments(p))

    @staticmethod
    def split(run, k):
        t = run.find(q('t')); s = t.text or ''
        if k <= 0 or k >= len(s): return
        right = copy.deepcopy(run)
        t.text = s[:k]; t.set(XML_SPACE, 'preserve')
        rt = right.find(q('t')); rt.text = s[k:]; rt.set(XML_SPACE, 'preserve')
        run.addnext(right)

    @staticmethod
    def top(run, p):
        """Outermost ancestor of run below the paragraph that is an ins wrapper, else the run itself."""
        node = run
        while node.getparent() is not p and node.getparent().tag == q('ins'):
            node = node.getparent()
        return node

    def ins_run(self, text, rpr_source):
        ins = self.mark('ins'); r = E.SubElement(ins, q('r'))
        rpr = rpr_source.find(q('rPr')) if rpr_source is not None else None
        if rpr is not None: r.append(copy.deepcopy(rpr))
        t = E.SubElement(r, q('t')); t.text = text; t.set(XML_SPACE, 'preserve')
        return ins

    # ------------------------------------------------------------ edits
    def edit(self, p, a, b, s):
        """Tracked replacement of characters [a, b) of the paragraph's visible text by s."""
        if b > a:
            for boundary in (b, a):
                for r, st, en in self.segments(p):
                    if st < boundary < en: self.split(r, boundary - st); break
            victims = [r for r, st, en in self.segments(p) if st >= a and en <= b and en > st]
            assert victims, (a, b, self.text(p)[a:b])
            last = None
            for r in victims:
                if r.getparent().tag == q('ins'):
                    raise ValueError('refusing to delete inside an inserted run')
                d = self.mark('del'); r.addprevious(d); d.append(r)
                r.find(q('t')).tag = q('delText'); last = d
            if s: last.addnext(self.ins_run(s, victims[0]))
            return
        if not s: return
        for r, st, en in self.segments(p):
            if st < a < en: self.split(r, a - st); break
        segs = [x for x in self.segments(p) if x[2] > x[1]]
        before = [r for r, st, en in segs if en == a]
        if before:
            self.top(before[-1], p).addnext(self.ins_run(s, before[-1]))
        elif segs:
            self.top(segs[0][0], p).addprevious(self.ins_run(s, segs[0][0]))
        else:
            ppr = p.find(q('pPr')); ins = self.ins_run(s, None)
            (ppr.addnext(ins) if ppr is not None else p.insert(0, ins))

    def replace(self, p, old, new, occurrence=0):
        text = self.text(p); start = -1
        for _ in range(occurrence + 1):
            start = text.find(old, start + 1)
            assert start >= 0, (old, text[:120])
        self.edit(p, start, start + len(old), new)

    def revise(self, p, new, gap=3):
        """Word-level diff between current text and new text, applied as tracked changes."""
        self.normalize(p)
        old = self.text(p)
        tok = lambda s: re.findall(r'\s+|\w+|[^\w\s]', s)
        a, b = tok(old), tok(new)
        ops = [o for o in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes()]
        blocks = []
        for tag, i1, i2, j1, j2 in ops:
            if tag == 'equal':
                continue
            if blocks:
                pi1, pi2, pj1, pj2 = blocks[-1]
                between = a[pi2:i1]
                if len(between) <= gap and ''.join(between).strip().__len__() <= 4:
                    blocks[-1] = (pi1, i2, pj1, j2); continue
            blocks.append((i1, i2, j1, j2))
        offs = [0]
        for t in a: offs.append(offs[-1] + len(t))
        for i1, i2, j1, j2 in reversed(blocks):
            self.edit(p, offs[i1], offs[i2], ''.join(b[j1:j2]))
        assert self.visible(p) == new, (self.visible(p), new)

    @staticmethod
    def visible(p):
        out = []
        for r in p.iter(q('r')):
            if any(x.tag == q('del') for x in r.iterancestors()): continue
            t = r.find(q('t'))
            if t is not None: out.append(t.text or '')
        return ''.join(out)

    # ------------------------------------------------------------ structure
    def mark_paragraph_inserted(self, p):
        ppr = p.find(q('pPr'))
        if ppr is None: ppr = E.Element(q('pPr')); p.insert(0, ppr)
        rpr = ppr.find(q('rPr'))
        if rpr is None:
            rpr = E.Element(q('rPr'))
            tail = [c for c in ppr if c.tag in (q('sectPr'), q('pPrChange'))]
            (tail[0].addprevious(rpr) if tail else ppr.append(rpr))
        rpr.insert(0, self.mark('ins'))

    def new_paragraph(self, after, text, template):
        p = E.Element(q('p'))
        ppr = template.find(q('pPr'))
        if ppr is not None:
            ppr = copy.deepcopy(ppr)
            for c in ppr.findall(q('sectPr')): ppr.remove(c)
            p.append(ppr)
        self.mark_paragraph_inserted(p)
        src = next((r for r in template.iter(q('r')) if r.find(q('t')) is not None), None)
        p.append(self.ins_run(text, src))
        after.addnext(p)
        return p

    def insert_row(self, tbl, after_index, texts):
        rows = tbl.findall(q('tr')); row = copy.deepcopy(rows[after_index])
        trpr = row.find(q('trPr'))
        if trpr is None:
            trpr = E.Element(q('trPr'))
            prex = row.find(q('tblPrEx'))
            (prex.addnext(trpr) if prex is not None else row.insert(0, trpr))
        trpr.append(self.mark('ins'))
        for tc, text in zip(row.findall(q('tc')), texts):
            paras = tc.findall(q('p'))
            for extra in paras[1:]: tc.remove(extra)
            p = paras[0]; src = next((r for r in p.iter(q('r')) if r.find(q('t')) is not None), None)
            src = copy.deepcopy(src) if src is not None else None
            for c in list(p):
                if c.tag != q('pPr'): p.remove(c)
            self.mark_paragraph_inserted(p)
            p.append(self.ins_run(text, src))
        rows[after_index].addnext(row)
        return row


def cell_paragraph(tbl, r, c):
    tc = tbl.findall(q('tr'))[r].findall(q('tc'))[c]
    ps = tc.findall(q('p'))
    return next((p for p in ps if Tracker.text(p)), ps[0])


def write_docx(source_docx, document_root, target):
    data = E.tostring(document_root, xml_declaration=True, encoding='UTF-8', standalone=True)
    with zipfile.ZipFile(source_docx) as zin, zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            zout.writestr(info, data if info.filename == 'word/document.xml' else zin.read(info.filename))
