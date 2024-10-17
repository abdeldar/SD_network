import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import time
from sklearn.inspection import permutation_importance
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from scipy import stats
import scipy.stats


def l_and_r_2_both(df_in):
    df_in["DNA"] = df_in["DNA_l"] + df_in["DNA_r"]
    df_in["LINE"] = df_in["LINE_l"] + df_in["LINE_r"]
    df_in["LTR"] = df_in["LTR_l"] + df_in["LTR_r"]
    df_in["SINE"] = df_in["SINE_l"] + df_in["SINE_r"]
    df_in["Low_complexity"] = df_in["Low_complexity_l"] + df_in["Low_complexity_r"]
    df_in["Retroposon"] = df_in["Retroposon_l"] + df_in["Retroposon_r"]
    df_in["Satellite"] = df_in["Satellite_l"] + df_in["Satellite_r"]
    df_in["Simple_repeat"] = df_in["Simple_repeat_l"] + df_in["Simple_repeat_r"]
    df_in["rRNA"] = df_in["rRNA_l"] + df_in["rRNA_r"]
    df_in["snRNA"] = df_in["snRNA_l"] + df_in["snRNA_r"]
    df_in["scRNA"] = df_in["scRNA_l"] + df_in["scRNA_r"]
    df_in["srpRNA"] = df_in["srpRNA_l"] + df_in["srpRNA_r"]
    df_in["tRNA"] = df_in["tRNA_l"] + df_in["tRNA_r"]
    df_in["RC"] = df_in["RC_l"] + df_in["RC_r"]  
    df_in["L1"] = df_in["L1_s_l"] + df_in["L1_s_r"]
    df_in["L2"] = df_in["L2_s_l"] + df_in["L2_s_r"]
    df_in["MIR"] = df_in["MIR_s_l"] + df_in["MIR_s_r"]
    df_in["Alu"] = df_in["Alu_s_l"] + df_in["Alu_s_r"]
    df_in["Real_satellite"] = df_in["Satellite_s_l"] + df_in["Satellite_s_r"]
    
    df_in["CG_frac"] = df_in[["CG_frac_l", "CG_frac_r"]].apply(lambda x: np.mean(x) if min(x)>=0 else max(x), axis=1)
    df_in["telocent_dist"] = df_in[["telo", "centro"]].apply(min, axis=1)
    
    df_in = df_in.drop(['DNA_l', 'LINE_l', 'LTR_l', 'SINE_l', 'Low_complexity_l', 'Retroposon_l', 'Satellite_l',
       'Simple_repeat_l', 'rRNA_l', 'snRNA_l', 'scRNA_l', 'srpRNA_l', 'tRNA_l',
       'RC_l', 'DNA_r', 'LINE_r', 'LTR_r', 'SINE_r', 'Low_complexity_r',
       'Retroposon_r', 'Satellite_r', 'Simple_repeat_r', 'rRNA_r', 'snRNA_r',
       'scRNA_r', 'srpRNA_r', 'tRNA_r', 'RC_r'], axis=1)
    df_in = df_in.drop(['L1_s_l', 'L2_s_l', 'MIR_s_l',
       'Alu_s_l', 'Satellite_s_l', 'L1_s_r', 'L2_s_r', 'MIR_s_r',
       'Alu_s_r', 'Satellite_s_r'], axis=1)
    df_in = df_in.drop(['CG_frac_l', 'CG_frac_r'], axis=1)
    df_in = df_in.drop(['used_coor_l_s', 'used_coor_l_e', 'used_coor_r_s', 'used_coor_r_e'], axis=1)
    df_in = df_in.drop(['repli_vari'], axis=1)
    return df_in


def plot_shifts_repeats(dfs, col_name, ax_in, sides, comp1="all", comp2="all", col_in="tab:blue"):
    arr = np.array([])
    for df in dfs:
        if comp1 == "noNs":
            df = df[df["gaps"] == 0]
        
        if comp2 == "small":
            df = df[df["jumps"] == 1]
        elif comp2 == "big":
            df = df[df["jumps"] > 3]
        arr = np.append(arr, sum(df[col_name])/(2*len(df[col_name])))
    ax_in.plot(range(-5, len(arr)-5), arr, color=col_in, linewidth=2, linestyle='-', marker='o', markersize=4)
    y_min, y_max = ax_in.get_ylim()
    rect = Rectangle((-15, -0.1), 15.5, y_max + 0.2, linewidth=0, facecolor='cornflowerblue', alpha=0.3, label='_nolegend_')
    ax_in.add_patch(rect)
    rect = Rectangle((0.5, -0.1), 30, y_max + 0.2, linewidth=0, facecolor='mediumaquamarine', alpha=0.3, label='_nolegend_')
    ax_in.add_patch(rect)
    rect = Rectangle((30.5, -0.1), 15, y_max + 0.2, linewidth=0, facecolor='lightgrey', alpha=0.9, label='_nolegend_')
    ax_in.add_patch(rect)
    ax_in.set_xlim(-5.2, 40.2)
    ax_in.set_ylim(y_min, y_max) 


def plot_several_repeats(pieces, repeat_names, titles, subplots_size, figure_size, sides, comp1="all", comp2="all"):
    fig, axs = plt.subplots(subplots_size[0], subplots_size[1], figsize=figure_size)
    fig.subplots_adjust(hspace=.3)
    fig.subplots_adjust(wspace=.25)
    for i, ax in enumerate(axs.flat):
        if i >= len(repeat_names):
            ax.axis('off')
        else:
            ax.set_yscale('log')
            ax.set_title(titles[i], fontsize=24)
            ax.set_xticks([-5, 0, 5, 10, 15, 20, 25, 30, 35, 40])
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=13, which="both")
            plot_shifts_repeats(pieces, repeat_names[i], ax, sides, comp1="noNs", comp2=comp2, col_in="tab:orange")
            plot_shifts_repeats(pieces, repeat_names[i], ax, sides, comp1=comp1, comp2=comp2, col_in="tab:blue")
            #handles, legends = ax.get_legend_handles_labels()
            #print(legends)
            #order = [0, 2, 1, 3, 4]
            ax.legend(['without gaps', 'all duplicated'], loc='lower center', fontsize=11, frameon=False)
    fig.show()


def significant_association_with_jumps(pie, perm=10, thr = 0.05):
    start = time.time()
    X = np.array(pie.loc[:, ~pie.columns.isin(['chr', 'coor_s', 'coor_e', 'ids', 'jumps'])], float)
    y = np.array(pie["jumps"], float)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    clf = RandomForestRegressor(max_depth=5, n_estimators=100, random_state=42, oob_score=True)
    clf.fit(X_train, y_train)
    imp_r = permutation_importance(clf, X_test, y_test)
    imps = np.array([imp_r.importances_mean])
    for i in range(perm):
        np.random.shuffle(y_test)
        np.random.shuffle(y_train)
        clf.fit(X_train, y_train)
        imp_p = permutation_importance(clf, X_test, y_test)
        imps = np.append(imps, [imp_p.importances_mean], axis=0)
    
    cols = np.array(pie.loc[:, ~pie.columns.isin(['chr', 'coor_s', 'coor_e', 'ids', 'jumps'])].columns)
    imps_arr = imps[0] - imps < 0
    imps_arr = imps_arr.astype(int)
    pvals = np.sum(imps_arr, axis=0)/(imps_arr.shape[0] - 1)
    pvals = stats.false_discovery_control(pvals)
    print(f"All columns: {cols}")
    print(f"Significant: {cols[(pvals < thr) & (imps[0] > 0)]}")
    print(f"Corresponding p-adjusted (BH): {pvals[(pvals < thr) & (imps[0] > 0)]}")


def breakpoints_enrichment(piece_out, pieces, mode="normal"):
    ms, stds, vals, labels, colors = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    cg_frac = piece_out.pop('CG_frac')
    piece_out.insert(23, 'CG_frac', cg_frac)
    cols = np.array(piece_out.columns[7:])
    cols = list(filter(lambda x: not x in ["intra_frac", "Retroposon", "Satellite", "rRNA", 
                                "snRNA", "scRNA", "srpRNA", "tRNA", "RC", "telocent_dist", "jumps"], cols))
    for colname in cols:
        if colname in ["genes", "cpgisl_in", "repli_in", "repli_deriv", "recomb_in", "dnase_in", "CG_frac_in", "CG_frac"]:
            norm_coef = 1
        else:
            norm_coef = 2
        if mode == "normal":
            val = np.mean(piece_out[colname])/norm_coef
            m_vs = [np.mean(pie.loc[(pie["CG_frac_in"] != 0.0), colname])/norm_coef for pie in pieces[-10:]]
        elif (mode == "nogaps") and (colname != "gaps"):
            val = np.mean(piece_out.loc[piece_out["gaps"] == 0.0, colname])/norm_coef
            m_vs = [np.mean(pie.loc[(pie["CG_frac_in"] != 0.0) & (pie["gaps"] == 0.0), colname])/norm_coef for pie in pieces[-10:]]
        elif colname == "gaps":
            ms = np.append(ms, 0.0)
            stds = np.append(stds, 0.0)
            vals = np.append(vals, 0.0)
            labels = np.append(labels, "gaps")
            colors = np.append(colors, 'k')
            continue
        m_v = np.mean(m_vs)
        std_v = np.std(m_vs)
        ci_small, ci1, ci2, ci_large = scipy.stats.norm.ppf([0.001/len(cols), 0.01/len(cols), 1 - 0.01/len(cols), 1 - 0.001/len(cols)], m_v, std_v) 
        if val > ci2:
            colors = np.append(colors, "r")
            if val > ci_large:
                print(f"{colname}, (sign. more): {val} {ci2}")
            else:
                print(f"{colname}, more: {val} {ci2}")
        elif val < ci1:
            colors = np.append(colors, "b")
            if val < ci_small:
                print(f"{colname}, (sign. less): {val} {ci1}")
            else:
                print(f"{colname}, less: {val} {ci1}")
        else:
            print(f"{colname}, NOTHING")
            continue
        ms = np.append(ms, m_v)
        stds = np.append(stds, std_v)
        vals = np.append(vals, val)
        labels = np.append(labels, colname)
    return [ms, stds, vals, labels, colors]
