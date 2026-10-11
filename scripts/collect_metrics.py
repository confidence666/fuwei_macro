#!/usr/bin/env python3
from pathlib import Path
import csv, re, sys

archive = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
benchmarks = ['bp','bp_be','ariane81','swerv_wrapper','swerv_wrapper43','mor1kx','or1200','ariane133','bp_multi','bp_multi57','isa_npu','dft68','ethernet','bp_fe','ariane136']
flows = ['global','global_incre']

def rel(p):
    try:
        return str(p.relative_to(archive))
    except ValueError:
        return str(p)

def first_error(text):
    for line in text.splitlines():
        if '[ERROR GPL-0305]' in line:
            continue
        if re.search(r'(\[ERROR|^Error:|Signal [0-9]+|segmentation fault|command failed|timeout)', line, re.I):
            return line.strip()[:240]
    return ''

def metric_after_marker(lines, marker, word):
    for i, line in enumerate(lines):
        if marker in line:
            for nxt in lines[i+1:i+20]:
                m = re.search(r'\b' + re.escape(word) + r'\s+(-?\d+(?:\.\d+)?)', nxt, re.I)
                if m:
                    return m.group(1)
    return ''

def last_metric(lines, word):
    val = ''
    pat = re.compile(r'\b' + re.escape(word) + r'\s+(-?\d+(?:\.\d+)?)', re.I)
    for line in lines:
        m = pat.search(line)
        if m:
            val = m.group(1)
    return val

def last_hpwl(lines):
    val = ''
    pat = re.compile(r'legalized HPWL\s+(-?\d+(?:\.\d+)?)', re.I)
    for line in lines:
        m = pat.search(line)
        if m:
            val = m.group(1)
    return val


rows = []
for design in benchmarks:
    for flow in flows:
        outdir = archive/'04_results'/design/flow
        log = outdir/f'run_{flow}.log'
        dfile = outdir/f'{design}_{flow}.def'
        odb = outdir/f'{design}_{flow}.odb'
        text = log.read_text(errors='replace') if log.exists() else ''
        lines = text.splitlines()
        status = 'pass'
        notes = ''
        if (outdir/'FAILED').exists():
            status = 'failed'
            notes = (outdir/'FAILED').read_text(errors='replace').strip()
            if dfile.exists() and odb.exists() and ('METRIC_' in text) and 'GPL-0305' in text and not first_error(text):
                status = 'pass_with_gpl0305'
                notes = 'completed with [ERROR GPL-0305] messages in log; final DEF/ODB and metrics exist'
        elif not log.exists():
            status = 'failed'; notes = 'missing log'
        elif first_error(text):
            status = 'failed'; notes = first_error(text)
        elif not dfile.exists() or not odb.exists():
            status = 'failed'; notes = 'missing expected DEF/ODB'
        if status == 'pass' and '[ERROR GPL-0305]' in text:
            status = 'pass_with_gpl0305'
            notes = 'completed with [ERROR GPL-0305] messages in log; final DEF/ODB and metrics exist'
        if flow == 'global_incre':
            wns = metric_after_marker(lines, 'METRIC_GLOBAL_INCRE_AFTER_WNS_BEGIN', 'wns') or metric_after_marker(lines, 'METRIC_AFTER_INCRE_WNS_BEGIN', 'wns') or last_metric(lines, 'wns')
            tns = metric_after_marker(lines, 'METRIC_GLOBAL_INCRE_AFTER_TNS_BEGIN', 'tns') or metric_after_marker(lines, 'METRIC_AFTER_INCRE_TNS_BEGIN', 'tns') or last_metric(lines, 'tns')
        else:
            wns = metric_after_marker(lines, 'METRIC_GLOBAL_WNS_BEGIN', 'wns') or last_metric(lines, 'wns')
            tns = metric_after_marker(lines, 'METRIC_GLOBAL_TNS_BEGIN', 'tns') or last_metric(lines, 'tns')
        hpwl = last_hpwl(lines)
        row = {
            'case': design,
            'flow': flow,
            'status': status,
            'wns': wns,
            'tns': tns,
            'hpwl': hpwl,
            'def': rel(dfile) if dfile.exists() else '',
            'odb': rel(odb) if odb.exists() else '',
            'log': rel(log) if log.exists() else '',
            'notes': notes,
        }
        rows.append(row)
        outdir.mkdir(parents=True, exist_ok=True)
        with (outdir/f'metrics_{flow}.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(row.keys()))
            writer.writeheader(); writer.writerow(row)

summary = archive/'04_results'/'summary_metrics.csv'
with summary.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['case','flow','status','wns','tns','hpwl','def','odb','log','notes'])
    writer.writeheader(); writer.writerows(rows)

passes = sum(1 for r in rows if r['status'].startswith('pass'))
print(f'wrote {summary} ({passes}/{len(rows)} flows completed)')
