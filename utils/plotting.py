"""
Plotting functions.
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotting_utils as plu
from circlify import _bubbles, circlify, Circle


##

        
def packed_circle_plot(
    df, covariate=None, ax=None, color='b', cmap=None, alpha=.5, linewidth=1.2,
    t_cov=.01, annotate=False, fontsize=6, ascending=False, fontcolor='white', 
    fontweight='normal'
    ):

    """
    Circle plot. Packed.
    """
    df = df.sort_values(covariate, ascending=False)
    circles = circlify(
        df[covariate].to_list(),
        show_enclosure=True, 
        target_enclosure=Circle(x=0, y=0, r=1)
    )
    lim = max(
        max(
            abs(c.x) + c.r,
            abs(c.y) + c.r,
        )
        for c in circles
    )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    
    if isinstance(color, str) and not color in df.columns:
        colors = { k : color for k in df.index }
    elif isinstance(color, str) and color in df.columns:
        c_cont = plu.create_palette(
            df.sort_values(color, ascending=True),
            color, cmap
        )
        colors = {}
        for name in df.index:
            colors[name] = c_cont[df.loc[name, color]]
    else:
        assert isinstance(color, dict)
        colors = color
        print('Try to use custom colors...')

    for name, circle in zip(df.index[::-1], circles): # Don't know why, but it reverses...
        x, y, r = circle
        ax.add_patch(
            plt.Circle((x, y), r*0.95, alpha=alpha, linewidth=linewidth, 
                fill=True, edgecolor=colors[name], facecolor=colors[name])
        )
        if annotate:
            cov = df.loc[name, covariate]
            if cov > t_cov:
                n = name if len(name)<=5 else name[:5]
                ax.annotate(
                    f'{n}: {df.loc[name, covariate]:.2f}', 
                    (x,y), 
                    va='center', ha='center', 
                    fontweight=fontweight, fontsize=fontsize, color=fontcolor, 
                )

    ax.axis('off')
    
    return ax


##