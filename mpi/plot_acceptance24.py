"""Plot the independently audited five complete MPI comparisons, not pilots."""
import argparse
import json
from pathlib import Path
import statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('evidence', type=Path)
    p.add_argument('audit', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    report = json.loads(a.audit.read_text())
    assert report['accepted'] and report['numerical_pass'] and report['nondestructive_snapshot_pass']
    candidate = list(report['timing_samples_seconds'].values())
    assert len(candidate) == 5
    base = []
    for i in range(1, 6):
        r = json.loads((a.evidence / 'runs' / ('r014-control%02d' % i) / 'result.json').read_text())
        assert r['passed'] and r['ranks'] == 32 and r['threads'] == 1 and r['steps'] == 4800
        base.append(r['elapsed'])
    fig, ax = plt.subplots(figsize=(8.5, 5.3), layout='constrained')
    colors = ['#426b9a', '#158363']
    meds = [statistics.median(base), statistics.median(candidate)]
    ax.bar([0, 1], meds, width=.54, color=colors, alpha=.18, edgecolor=colors, linewidth=1.5)
    for i, (x, y) in enumerate(zip(base, candidate)):
        dx = (i - 2) * .055
        ax.plot([dx, 1+dx], [x, y], color='#ccd2d6', linewidth=.8, zorder=2)
    for j, values in enumerate((base, candidate)):
        ax.scatter([j+(i-2)*.055 for i in range(5)], values, s=30, color=colors[j], zorder=4)
        ax.hlines(meds[j], j-.23, j+.23, colors=colors[j], linewidth=2.5, zorder=3)
        ax.text(j, max(values)+.85, f'median {meds[j]:.3f} s', ha='center', color=colors[j], weight='bold')
    ax.axhline(report['target_seconds'], color='#b34738', linestyle='--', linewidth=1.15)
    ax.text(-.43, report['target_seconds']+.12, '24x target: 13.565 s',
            color='#b34738', fontsize=9, ha='left', bbox=dict(facecolor='white', edgecolor='none', pad=2))
    ax.set(xticks=[0, 1], xticklabels=['Previous MPI (r004b)', 'Accepted MPI (r014)'],
           ylabel='Complete run wall time (seconds; lower is better)', xlim=(-.5, 1.5), ylim=(0, 20.3))
    ax.set_title('sECOM on p3: 32 MPI ranks x 1 thread\nFive new interleaved comparisons', loc='left', weight='bold')
    ax.grid(axis='y', alpha=.18)
    ax.spines[['top', 'right']].set_visible(False)
    fig.supxlabel(f"{report['frozen_median_speedup']:.2f}x vs original frozen serial; "
                  f"{report['same_day_median_speedup']:.2f}x vs same-day original serial\n"
                  'Startup and required output included. Not same-revision 1-MPI scaling.', fontsize=10)
    assert not a.output.exists()
    fig.savefig(a.output, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
