import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"

if __name__ == "__main__":
    """
    Social Welfare Maximization Plot 1
    """
    N = 1000
    UL = 1
    x = np.linspace(0, UL, N)
    y1 = x
    y2 = -x + np.ones(N)

    fig, axes = plt.subplots(1, 2, figsize=(
        14, 6), facecolor='w', sharey=True, sharex=True)

    axes[0].plot(x, y2, label='Demand')
    axes[0].plot(x, y1, label='Supply')

    axes[0].plot([0.5], [0.5], 'ro', markersize=10)
    axes[0].axhline(y=0.5, xmax=0.5, linestyle='--', alpha=0.6, color='red')
    axes[0].fill_between(x=x, y1=0.5, y2=y2, where=x < 0.5, alpha=0.2)
    axes[0].fill_between(
        x=x,
        y1=y1,
        y2=0.5,
        where=x < 0.5,
        alpha=0.2,
        color='tab:orange')

    axes[0].set_ylabel('Normalized Price', fontsize=14)
    axes[0].set_xlabel('Normalized Quantity', fontsize=14)
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].text(0.05, 0.6, "Consumer Surplus", fontsize=12)
    axes[0].text(0.05, 0.4, "Producer Surplus", fontsize=12)

    axes[0].text(0.8, 0.8, 'Supply',
                 rotation=45, fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')
    axes[0].text(0.8, 0.2, 'Demand',
                 rotation=-40, fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')

    p = 0.3
    axes[1].plot(x, y2, label='Demand')
    axes[1].plot(x, y1, label='Supply')

    axes[1].plot([p], [p], 'ro', markersize=10)
    axes[1].axhline(y=p, xmax=p, linestyle='--', alpha=0.6, color='red')
    axes[1].fill_between(x=x, y1=p, y2=y2, where=x < p, alpha=0.2)
    axes[1].fill_between(
        x=x,
        y1=y1,
        y2=p,
        where=x < p,
        alpha=0.2,
        color='tab:orange')
    axes[1].fill_between(
        x=x,
        y1=y1,
        y2=y2,
        where=(
            (x > p) & (
                x < 0.5)),
        alpha=0.2,
        color='gray')

    axes[1].set_xlabel('Normalized Quantity', fontsize=14)
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1)
    axes[1].text(0.05, 0.5, "Consumer \nSurplus", fontsize=12)
    axes[1].text(0.05, 0.2, "Producer \nSurplus", fontsize=12)
    axes[1].text(0.32, 0.45, "Lost \nSurplus", fontsize=12)

    axes[1].text(0.8, 0.8, 'Supply',
                 rotation=45, fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')
    axes[1].text(0.8, 0.2, 'Demand',
                 rotation=-40, fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')

    axes[0].text(0.01, 0.95, "a)", fontsize=14, bbox=dict({'facecolor': 'w'}))
    axes[1].text(0.01, 0.95, "b)", fontsize=14, bbox=dict({'facecolor': 'w'}))
    fig.suptitle('Social Welfare Maximization', fontsize=16)
    plt.tight_layout()
    plt.savefig('../docs/figures/social_max.pgf')

    """
    Social Welfare Maximization Plot 2: Elasticity
    """
    N = 1000
    UL = 1
    x = np.linspace(0, UL, N)
    y1 = x
    y2 = -x + np.ones(N)
    a1 = b1 = 0.5
    a2 = b2 = 0.0
    md = -1000
    ms = 0.0
    def d(x): return md * (x - a1) + b1
    def s(x): return ms * (x - a1) + b1
    y3 = d(x)
    y4 = s(x)

    fig, axes = plt.subplots(1, 2, figsize=(
        14, 6), facecolor='w', sharey=True, sharex=True)

    axes[0].plot(x, y2, label='Demand')
    axes[0].plot(x, y1, label='Supply')
    axes[0].axhline(y=0.5, xmax=0.5, linestyle='--', alpha=0.6, color='red')

    axes[0].plot([0.5], [0.5], 'ro', markersize=10)
    axes[0].fill_between(x=x, y1=0.5, y2=y2, where=x < 0.5, alpha=0.2)
    axes[0].fill_between(
        x=x,
        y1=y1,
        y2=0.5,
        where=x < 0.5,
        alpha=0.2,
        color='tab:orange')

    axes[0].set_ylabel('Normalized Price', fontsize=14)
    axes[0].set_xlabel('Normalized Quantity', fontsize=14)
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].text(0.05, 0.6, "Consumer Surplus", fontsize=12)
    axes[0].text(0.05, 0.4, "Producer Surplus", fontsize=12)

    axes[0].text(0.8, 0.8, 'Supply',
                 rotation=45, fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')
    axes[0].text(0.8, 0.2, 'Demand',
                 rotation=-40, fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')

    p = 0.5
    axes[1].plot(x, y3, label='Demand')
    axes[1].plot(x, y4, label='Supply')

    idx = np.argwhere(np.diff(np.sign(y4 - y3))).flatten()
    eqp = 0.5
    p = y4[idx]
    q = x[idx]

    axes[1].plot(q, p, 'ro', markersize=10)
    axes[1].axhline(y=p, xmax=q[0], linestyle='--', alpha=0.6, color='red')

    axes[1].fill_between(x=x, y1=p, y2=y3, where=x < q, alpha=0.2)
    axes[1].fill_between(
        x=x,
        y1=0,
        y2=p,
        where=x < q,
        alpha=0.2,
        color='tab:orange')

    axes[1].set_xlabel('Normalized Quantity', fontsize=14)
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1)
    axes[1].text(0.05, 0.6, "Consumer Surplus", fontsize=12)
    axes[1].text(0.05, 0.4, "Producer Surplus", fontsize=12)

    axes[1].text(0.8, p - 0.02, 'Supply',
                 rotation=0, fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')
    axes[1].text(0.47, 0.2, 'Demand',
                 rotation=-90, fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')

    axes[0].text(0.01, 0.95, "a)", fontsize=14, bbox=dict({'facecolor': 'w'}))
    axes[1].text(0.01, 0.95, "b)", fontsize=14, bbox=dict({'facecolor': 'w'}))
    fig.suptitle("Price Elasticity", fontsize=16)
    plt.tight_layout()
    plt.savefig('../docs/figures/elasticity.pgf')
