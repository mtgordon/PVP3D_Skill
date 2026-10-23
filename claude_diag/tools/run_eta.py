"""Run progress + ETA (2026-10-23): t, gain since the last call, rate, estimated t=1 and t=1.5 times (hold at 0.9 x the rate).
usage: py -3.10 run_eta.py RUN [RUN ...]   (state in runs/_eta_state.json)"""
import re, json, os, datetime, sys
RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
st = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs', '_eta_state.json')
names = sys.argv[1:]
prev = json.load(open(st)) if os.path.exists(st) else {}
now = datetime.datetime.now(); cur = {}
print(now.strftime('%H:%M'))
for r in names:
    log = open(os.path.join(RUNS, r, r + '.log'), errors='ignore').read()
    ts = re.findall(r'converged at time : ([\d.e+-]+)', log); t = float(ts[-1]) if ts else 0.0
    done = 'N O R M A L' in log; cur[r] = (t, now.isoformat())
    msg = f'{r[:5]:6s} t {t:.3f}'
    if done:
        print(msg + '  NORMAL TERMINATION'); continue
    if r in prev:
        t0, w0 = prev[r]; dt = (now - datetime.datetime.fromisoformat(w0)).total_seconds() / 3600
        rate = (t - t0) / dt if dt > 0 else 0
        if rate > 1e-4:
            hold = rate * 0.9
            a = (max(0, 1 - t) / rate); b = (1.5 - max(t, 1)) / hold
            t1 = (now + datetime.timedelta(hours=a)).strftime('%H:%M') if t < 1 else 'passed'
            msg += f'  +{t - t0:.3f} in {dt * 60:.0f} min ({rate:.3f} t/h)  t=1 {t1}  t=1.5 ~{(now + datetime.timedelta(hours=a + b)).strftime("%H:%M")}'
        else:
            msg += f'  NO PROGRESS in {dt * 60:.0f} min'
    print(msg)
json.dump(cur, open(st, 'w'))
