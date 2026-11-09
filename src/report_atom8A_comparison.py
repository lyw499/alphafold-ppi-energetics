"""Raw and conservative-QC distributions, with no count-based outlier filtering."""
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/atom8A-matplotlib')
import numpy as np
import pandas as pd

from atom_contact_common import CUTOFF, DEFINITION

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/contact_comparison_atom8A'
MAIN = [('Pedro', 'AF3'), ('SAFES', 'experimental'), ('Mgen_reference', 'experimental'), ('FoldBench', 'experimental')]
LABELS = {('Pedro', 'AF3'): 'Pedro (AF3)', ('SAFES', 'experimental'): 'SAFES (experimental)',
          ('Mgen_reference', 'experimental'): 'M. genitalium candidates (experimental)',
          ('FoldBench', 'experimental'): 'FoldBench (experimental)',
          ('SAFES', 'AF3'): 'SAFES (AF3)', ('FoldBench', 'AF3'): 'FoldBench (AF3)'}


def stats(values):
    x = pd.Series(values).dropna()
    return dict(mean_contact_count=x.mean(), median_contact_count=x.median(), std_contact_count=x.std(ddof=1),
                q1=x.quantile(.25), q3=x.quantile(.75), min=x.min(), max=x.max())


def ecdf_distance(a, b):
    a, b = np.sort(a), np.sort(b)
    if not len(a) or not len(b):
        return np.nan
    grid = np.unique(np.concatenate([a, b]))
    return float(np.max(np.abs(np.searchsorted(a, grid, side='right') / len(a) - np.searchsorted(b, grid, side='right') / len(b))))


def plot_ecdf(ax, values, label, color=None, linestyle='-'):
    x = np.sort(np.asarray(values, dtype=float))
    if len(x):
        ax.step(x, np.arange(1, len(x) + 1) / len(x), where='post', label=f'{label} (n={len(x):,})', color=color, linestyle=linestyle)


def load_data(allow_partial=False):
    manifest = json.loads((BASE / 'pedro_manifest.json').read_text())
    complete = [s for s in manifest if (BASE / 'pedro_parts' / s['source_file']).with_suffix('.done.json').exists()]
    if len(complete) != len(manifest) and not allow_partial:
        raise RuntimeError(f'Pedro incomplete: {len(complete)}/{len(manifest)} archives. No final comparison written.')
    frames = [pd.read_csv(BASE / 'pedro_parts' / s['source_file']) for s in complete]
    for name in ('safes_experimental', 'safes_af3', 'mgen_reference_experimental', 'foldbench_experimental', 'foldbench_af3'):
        frames.append(pd.read_csv(BASE / (name + '.csv')))
    data = pd.concat(frames, ignore_index=True)
    assert data.contact_definition.eq(DEFINITION).all()
    assert data.cutoff_angstrom.eq(CUTOFF).all()
    if 'malformed_chain_ids' in data:
        assert not (data.malformed_duplicate_identity.fillna(False) & data.malformed_chain_ids.isna()).any(), 'Unverified parser-warning QC needs review'
    for col in ('raw_included', 'qc_included', 'included_main_distribution', 'included_qc_nonzero_distribution'):
        assert data[col].isin([True, False]).all(), col
    assert (data.loc[data.raw_included, 'contact_count'] >= 0).all()
    return data, dict(pedro_archives_complete=len(complete), pedro_archives_expected=len(manifest), final=len(complete) == len(manifest))


def pair_level(pedro, out):
    summaries, versions = [], {}
    for version, mask in [('raw', pedro.raw_included), ('qc_filtered', pedro.qc_included)]:
        subset = pedro[mask]
        jobs = subset.groupby(['job_id', 'protein_pair'], as_index=False).agg(
            contact_count=('contact_count', 'mean'), n_archived_observations=('contact_count', 'size'))
        pairs = jobs.groupby('protein_pair', as_index=False).agg(
            contact_count=('contact_count', 'mean'), n_jobs=('job_id', 'size'),
            n_archived_observations=('n_archived_observations', 'sum'))
        versions[version] = pairs
        pairs.to_csv(out / ('pedro_pair_level_' + version + '.csv'), index=False)
        summaries.append(dict(version=version, observation_unit='protein_pair_two_stage_mean', n_total=len(pairs),
                              n_zero_contacts=int((pairs.contact_count == 0).sum()), **stats(pairs.contact_count)))
    pd.DataFrame(summaries).to_csv(out / 'pedro_pair_level_summary.csv', index=False)
    # Legacy pair-level exclusion is a separate sensitivity, not the default QC.
    legacy = pd.read_csv(ROOT / 'results/pair_level_contacts_8A.csv', usecols=['protein_pair', 'has_severe_overlap'])
    bad = set(legacy.loc[legacy.has_severe_overlap, 'protein_pair'])
    alternate = versions['raw'][~versions['raw'].protein_pair.isin(bad)]
    alternate.to_csv(out / 'pedro_pair_level_exclude_any_legacy_overlap.csv', index=False)
    return versions


def make_report(data, status, out):
    out.mkdir(parents=True, exist_ok=True)
    (out / 'figures').mkdir(exist_ok=True)
    data.to_csv(out / 'comparison_long.csv', index=False)
    data[data.raw_included].to_csv(out / 'raw_minimally_filtered.csv', index=False)
    data[data.qc_included].to_csv(out / 'qc_filtered.csv', index=False)
    data[~data.raw_included | (data.raw_included & ~data.qc_included)].to_csv(out / 'exclusion_audit.csv', index=False)
    summary, changes, flag_rows, nonzero_summary = [], [], [], []
    groups = {(a, b): d.copy() for (a, b), d in data.groupby(['dataset', 'structure_type'], sort=False)}
    for (dataset, kind), group in groups.items():
        raw, qc = group[group.raw_included], group[group.qc_included]
        group.to_csv(out / f'{dataset}_{kind}_all_rows.csv', index=False)
        for version, subset in [('raw', raw), ('qc_filtered', qc)]:
            x = subset.contact_count
            subset.to_csv(out / f'{dataset}_{kind}_{version}.csv', index=False)
            row = dict(dataset=dataset, structure_type=kind, version=version,
                       observation_unit='structure_chain_pair', n_total=len(group),
                       n_structures_total=group.structure_id.nunique(),
                       n_structures_calculable=raw.structure_id.nunique(),
                       n_not_calculable=len(group) - len(raw), n_calculable=len(raw),
                       n_qc_excluded=len(raw) - len(qc), n_included=len(subset),
                       n_zero_contacts=int((x == 0).sum()),
                       pct_zero_contacts=100 * (x == 0).mean() if len(x) else np.nan,
                       n_nonzero=int((x > 0).sum()), contact_definition=DEFINITION, cutoff_angstrom=CUTOFF,
                       **stats(x))
            summary.append(row)
            nonzero_summary.append(dict(dataset=dataset, structure_type=kind, version=version,
                                        observation_scope='nonzero_display_only', n_included=int((x > 0).sum()),
                                        **stats(x[x > 0])))
        nonzero_raw, nonzero_qc = raw.loc[raw.contact_count > 0, 'contact_count'], qc.loc[qc.contact_count > 0, 'contact_count']
        changes.append(dict(dataset=dataset, structure_type=kind, n_raw=len(raw), n_qc=len(qc),
                            n_removed=len(raw) - len(qc), percent_removed=100 * (len(raw) - len(qc)) / len(raw) if len(raw) else np.nan,
                            median_raw=raw.contact_count.median(), median_qc=qc.contact_count.median(),
                            median_shift=qc.contact_count.median() - raw.contact_count.median(),
                            mean_raw=raw.contact_count.mean(), mean_qc=qc.contact_count.mean(),
                            nonzero_median_raw=nonzero_raw.median(), nonzero_median_qc=nonzero_qc.median(),
                            max_raw=raw.contact_count.max(), max_qc=qc.contact_count.max(),
                            maximum_ecdf_distance_all=ecdf_distance(raw.contact_count, qc.contact_count),
                            maximum_ecdf_distance_nonzero=ecdf_distance(nonzero_raw, nonzero_qc)))
        for flag in sorted({f for text in group.qc_flag.fillna('') for f in text.split(';') if f}):
            n = group.qc_flag.fillna('').str.split(';').apply(lambda fs: flag in fs).sum()
            flag_rows.append(dict(dataset=dataset, structure_type=kind, qc_flag=flag, n_rows=int(n),
                                  used_for_qc_exclusion=flag in ('severe_overlap_lt_0.5A', 'malformed_duplicate_identity')))
    summary = pd.DataFrame(summary)
    changes = pd.DataFrame(changes)
    summary.to_csv(out / 'comparison_summary.csv', index=False)
    pd.DataFrame(nonzero_summary).to_csv(out / 'nonzero_display_summary.csv', index=False)
    summary[summary.dataset != 'Pedro'].to_csv(out / 'completed_control_summary.csv', index=False)
    changes.to_csv(out / 'qc_effect_comparison.csv', index=False)
    pd.DataFrame(flag_rows).to_csv(out / 'qc_flag_counts.csv', index=False)
    ref = groups[('Mgen_reference', 'experimental')]
    taxonomy_rows = []
    if 'reference_taxonomy_scope' in ref:
        for scope, subset in [('all_candidate_chain_pairs', ref),
                              ('at_least_one_Mgen_chain', ref[ref.reference_taxonomy_scope.isin(['both_Mgen', 'one_Mgen_partner'])]),
                              ('both_chains_Mgen', ref[ref.reference_taxonomy_scope == 'both_Mgen']),
                              ('neither_chain_Mgen', ref[ref.reference_taxonomy_scope == 'neither_Mgen'])]:
            subset = subset[subset.raw_included]
            taxonomy_rows.append(dict(scope=scope, n_interfaces=len(subset), n_zero_contacts=int((subset.contact_count == 0).sum()), **stats(subset.contact_count)))
        pd.DataFrame(taxonomy_rows).to_csv(out / 'reference_taxonomy_sensitivity.csv', index=False)
    if ('Pedro', 'AF3') in groups:
        pair_versions = pair_level(groups[('Pedro', 'AF3')], out)
    else:
        pair_versions = {}
    # Metadata-based sensitivity only: do not recompute structures or change raw.
    fb = groups[('FoldBench', 'AF3')]
    stable = fb[fb.raw_included & ~fb.qc_flag.fillna('').str.contains('chain_mapping_tie_changes_atom_count')]
    stable.to_csv(out / 'foldbench_af3_mapping_stable_sensitivity.csv', index=False)
    pd.DataFrame([dict(scope='raw_AF3', n=len(fb[fb.raw_included]), **stats(fb.loc[fb.raw_included, 'contact_count'])),
                  dict(scope='exclude_tie_sensitive_assignment_only', n=len(stable), **stats(stable.contact_count))]).to_csv(out / 'foldbench_mapping_sensitivity_summary.csv', index=False)
    plot_all(groups, pair_versions, out, status)
    (out / 'completion_status.json').write_text(json.dumps(status, indent=2) + '\n')
    write_text_report(summary, changes, status, out)
    validate(data, summary, out)
    print(summary[['dataset', 'structure_type', 'version', 'n_included', 'n_zero_contacts', 'median_contact_count', 'mean_contact_count']].to_string(index=False), flush=True)


def save(fig, out, name):
    import matplotlib.pyplot as plt
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out / 'figures' / (name + '.' + ext), dpi=180)
    plt.close(fig)


def plot_all(groups, pair_versions, out, status):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    prefix = '' if status['final'] else 'PARTIAL PEDRO — '
    for key, group in groups.items():
        raw = group.loc[group.included_main_distribution, 'contact_count']
        qc = group.loc[group.included_qc_nonzero_distribution, 'contact_count']
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        plot_ecdf(ax, raw, 'Raw/minimally filtered', '#2878B5')
        plot_ecdf(ax, qc, 'Conservative QC', '#D66B27', '--')
        ax.set_xscale('log')
        zero = int(((group.contact_count == 0) & group.raw_included).sum())
        individual_prefix = prefix if key[0] == 'Pedro' else ''
        ax.set(title=individual_prefix + LABELS[key] + f'\nNonzero interfaces only; raw zeros retained in tables: {zero:,}',
               xlabel='Heavy-atom pair contacts within 8 Å (log scale)', ylabel='Cumulative fraction', ylim=(0, 1.02))
        ax.legend(frameon=False)
        save(fig, out, key[0] + '_' + key[1] + '_raw_vs_qc_nonzero')
    colors = ['#2878B5', '#D66B27', '#27835D', '#8764A3']
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharex=True, sharey=True)
    for ax, version, column in zip(axes, ['Raw/minimally filtered', 'Conservative QC'], ['included_main_distribution', 'included_qc_nonzero_distribution']):
        for key, color in zip(MAIN, colors):
            if key in groups:
                values = groups[key].loc[groups[key][column], 'contact_count']
                plot_ecdf(ax, values, LABELS[key], color)
        ax.set_xscale('log')
        ax.set(title=version, xlabel='Heavy-atom pair contacts within 8 Å (log scale)', ylabel='Cumulative fraction', ylim=(0, 1.02))
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle(prefix + 'Four datasets: nonzero interfaces only; full contact range, zero counts preserved in tables')
    save(fig, out, 'four_datasets_raw_and_qc_nonzero_ecdf')
    # Companion figure includes zeros in both the display and population.
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharex=True, sharey=True)
    for ax, version, column in zip(axes, ['Raw/minimally filtered', 'Conservative QC'], ['raw_included', 'qc_included']):
        for key, color in zip(MAIN, colors):
            if key in groups:
                values = np.log10(groups[key].loc[groups[key][column], 'contact_count'] + 1)
                plot_ecdf(ax, values, LABELS[key], color)
        ax.set(title=version, xlabel='log10(heavy-atom pair contact count + 1)', ylabel='Cumulative fraction', ylim=(0, 1.02))
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle(prefix + 'Four datasets: all calculable interfaces, INCLUDING zeros')
    save(fig, out, 'four_datasets_raw_and_qc_all_including_zeros_ecdf')
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
    for ax, title, column in zip(axes, ['Raw/minimally filtered', 'Conservative QC'], ['included_main_distribution', 'included_qc_nonzero_distribution']):
        keys = [key for key in MAIN if key in groups]
        values = [groups[key].loc[groups[key][column], 'contact_count'].to_numpy() for key in keys]
        ax.boxplot(values, tick_labels=[key[0].replace('_', '\n') for key in keys], showfliers=True,
                   flierprops={'markersize': 2, 'alpha': .25})
        ax.set_yscale('log')
        ax.set(title=title, ylabel='8 Å heavy-atom pair contacts (log scale)')
    fig.suptitle(prefix + 'Nonzero distributions; outliers shown, no range trimming')
    save(fig, out, 'four_datasets_raw_and_qc_nonzero_boxplots')
    if pair_versions:
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        for version, df in pair_versions.items():
            plot_ecdf(ax, df.loc[df.contact_count > 0, 'contact_count'], version)
        ax.set_xscale('log')
        ax.set(title='Pedro protein-pair sensitivity: two-stage averaging\nQC applied to observations BEFORE averaging; nonzero only',
               xlabel='Mean heavy-atom pair contact count (log scale)', ylabel='Cumulative fraction')
        ax.legend(frameon=False)
        save(fig, out, 'pedro_pair_level_raw_vs_qc_nonzero')
    for dataset in ('SAFES', 'FoldBench'):
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        for kind in ('experimental', 'AF3'):
            g = groups[(dataset, kind)]
            plot_ecdf(ax, g.loc[g.included_main_distribution, 'contact_count'], kind)
        ax.set_xscale('log')
        ax.set(title=dataset + ': experimental and AF3 kept separate\nRaw nonzero chain-pair distributions; populations may differ',
               xlabel='8 Å heavy-atom pair contacts (log scale)', ylabel='Cumulative fraction')
        ax.legend(frameon=False)
        save(fig, out, dataset + '_experimental_vs_AF3_raw_nonzero')


def validate(data, summary, out):
    valid = data[data.raw_included]
    assert valid.contact_count.notna().all()
    assert np.all(valid.contact_count == valid.contact_count.astype('int64'))
    assert (data.qc_included <= data.raw_included).all()
    assert data.included_main_distribution.equals(data.raw_included & (data.contact_count > 0))
    assert data.included_qc_nonzero_distribution.equals(data.qc_included & (data.contact_count > 0))
    assert ((summary.n_included == summary.n_nonzero + summary.n_zero_contacts)).all()
    assert ((summary.n_total == summary.n_calculable + summary.n_not_calculable)).all()
    report = dict(status='passed', n_rows=len(data), n_calculable=len(valid),
                  checks=['one atom-pair definition in every dataset', 'integer structure-level counts',
                          'raw and QC retained separately', 'zero counts included in statistics',
                          'nonzero figure flags explicitly exclude only zeros from each version',
                          'summary population arithmetic'])
    (out / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')


def write_text_report(summary, changes, status, out):
    lines = ['# 四数据集 8 Å 重原子对计数比较', '',
             ('完整结果。' if status['final'] else '**Pedro 尚未完成；本目录为调试预览，不作为最终比较。**'), '',
             '原始/最小过滤结果与 QC 结果并列；QC 不自动优先。统计包含零接触，主非零图明确排除零值，另有包含零值的 ECDF。', '',
             '## 统计（包含零值）', '',
             '| 数据集 | 类型 | 版本 | 可用接口 | 零接触 | 零值% | 均值 | 中位数 | Q1 | Q3 | 最大值 |',
             '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, r in summary.iterrows():
        lines.append(f'| {r.dataset} | {r.structure_type} | {r.version} | {r.n_included} | {r.n_zero_contacts} | {r.pct_zero_contacts:.2f} | {r.mean_contact_count:.2f} | {r.median_contact_count:.2f} | {r.q1:.2f} | {r.q3:.2f} | {r["max"]:.2f} |')
    lines += ['', '## QC 的实际影响', '']
    for _, r in changes.iterrows():
        if r.n_removed == 0:
            lines.append(f'- {r.dataset} / {r.structure_type}：没有案例满足当前明确的 QC 排除规则；raw 与 QC 分布完全相同。')
        else:
            lines.append(f'- {r.dataset} / {r.structure_type}：排除 {r.n_removed} 个接口（{r.percent_removed:.3f}%）；包含零值的中位数 {r.median_raw:.2f} → {r.median_qc:.2f}，非零中位数 {r.nonzero_median_raw:.2f} → {r.nonzero_median_qc:.2f}，均值 {r.mean_raw:.2f} → {r.mean_qc:.2f}，最大值 {r.max_raw:.2f} → {r.max_qc:.2f}。完整/非零 ECDF 最大纵向差分别为 {r.maximum_ecdf_distance_all:.5f} / {r.maximum_ecdf_distance_nonzero:.5f}。不使用人为“明显变化”阈值，按这些实际变化判断主体或尾部是否受影响。')
    lines += ['', '## 解释限制', '',
              '- Pedro 主比较的单位是归档结构中的链对观测，重跑也保留；另提供按旧方法先平均重跑、再平均任务的蛋白对级结果。QC 先移除问题结构观测，再重新聚合。',
              '- SAFES 和参考 assembly 没有指定结合伙伴表，枚举全部无序蛋白/肽链对；FoldBench 使用官方接口清单。候选链对的选择机制不同。',
              '- 主四组选择 SAFES 实验和 FoldBench 实验作为参考，AF3 版本单独保留。SAFES 多链亚基界面不等同于指定结合伙伴界面。',
              '- 同一结构中的不同链对不是独立实验；重复结构/assembly 的代表性和抽样偏差未通过 QC 人为修正。',
              '- 参考候选包含不同物种来源的蛋白。来源范围属于数据集定义敏感性，不是结构质量过滤；见 reference_taxonomy_sensitivity.csv。物种标注使用基因来源，不把表达宿主当作蛋白来源。',
              '- 不因接触数大、低 ipTM、全局 RMSD 大、长度差异或零接触而过滤。链映射不确定性保留标记并另给 FoldBench 敏感性版本。',
              '- 解析/原子规则及所有脚本和文件说明见 METHODS.md；逐案例原因见 exclusion_audit.csv 与 qc_flag_counts.csv。']
    (out / 'REPORT.md').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--allow-partial', action='store_true', help='Write only to preview/; never overwrite final comparison')
    args = parser.parse_args()
    try:
        data, status = load_data(args.allow_partial)
    except Exception as exc:
        if not args.allow_partial:
            (BASE / 'completion_status.json').write_text(json.dumps(dict(final=False, status='blocked_or_failed', error=str(exc)), indent=2) + '\n')
        raise
    target = BASE / 'preview' if args.allow_partial else BASE
    make_report(data, status, target)
