# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# **************************************************************************

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from pyworkflow.protocol import params
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER

from ..protocols import ProtSpatialPPI


class ProtSpatialPPIViewer(ProtocolViewer):
    _label = 'viewer'
    _targets = [ProtSpatialPPI]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

    # =========================================================================
    # Parameters
    # =========================================================================

    def _defineParams(self, form):
        form.addSection(label='Interaction score visualization')

        form.addParam(
            'threshold',
            params.FloatParam,
            label='Score threshold: ',
            default=0.0,
            help='Interaction probability threshold (0-1).'
        )

        form.addParam(
            'filterMode',
            params.EnumParam,
            choices=['Above threshold (≥)', 'Below threshold (≤)'],
            default=0,
            display=params.EnumParam.DISPLAY_HLIST,
            label='Display combinations: ',
            help='Select which combinations to display according to their '
                 'interaction probability.'
        )

        form.addParam(
            'showHeatmap',
            params.LabelParam,
            label='Display interaction probability heatmap: ',
            help='Display interaction probabilities for all selected '
                 'structure or sequence combinations.'
        )

    # =========================================================================
    # Visualization
    # =========================================================================

    def _getVisualizeDict(self):
        return {
            'showHeatmap': self._showHeatmap
        }

    def _validate(self):
        return []

    # =========================================================================
    # Results
    # =========================================================================

    def _getResults(self):
        resultsFile = self.protocol._getExtraPath('spatialppi_results.csv')

        if not os.path.exists(resultsFile):
            raise FileNotFoundError(
                f'SpatialPPI results file not found: {resultsFile}'
            )

        df = pd.read_csv(resultsFile)
        df.columns = df.columns.str.strip()

        scoreColumn = 'Interaction probability'

        if scoreColumn not in df.columns:
            raise ValueError(
                f'Could not find {scoreColumn}. '
                f'Available columns: {list(df.columns)}'
            )

        # Identify input type from CSV columns
        if {'Structure 1', 'Chain 1', 'Structure 2', 'Chain 2'}.issubset(df.columns):

            firstColumn = 'Protein 1'
            secondColumn = 'Protein 2'

            # Combine structure names and chain identifiers
            df[firstColumn] = (
                    df['Structure 1'].astype(str) + '_' +
                    df['Chain 1'].astype(str)
            )

            df[secondColumn] = (
                    df['Structure 2'].astype(str) + '_' +
                    df['Chain 2'].astype(str)
            )

        elif {'Protein 1', 'Protein 2'}.issubset(df.columns):

            firstColumn = 'Protein 1'
            secondColumn = 'Protein 2'

        else:
            raise ValueError(
                'Could not identify the interaction partner columns. '
                f'Available columns: {list(df.columns)}'
            )

        # Convert probabilities to numeric values
        df[scoreColumn] = pd.to_numeric(
            df[scoreColumn], errors='coerce'
        )

        df = df.dropna(
            subset=[firstColumn, secondColumn, scoreColumn]
        )

        return df, firstColumn, secondColumn, scoreColumn

    # =========================================================================
    # Heatmap
    # =========================================================================

    def _showHeatmap(self, paramName=None):

        df, firstColumn, secondColumn, scoreColumn = self._getResults()

        threshold = self.threshold.get()

        # Apply threshold
        if self.filterMode.get() == 0:
            df = df[df[scoreColumn] >= threshold]
            filterLabel = f'Probability ≥ {threshold:.3f}'
        else:
            df = df[df[scoreColumn] <= threshold]
            filterLabel = f'Probability ≤ {threshold:.3f}'

        if df.empty:
            print('No interaction combinations satisfy the selected '
                  'threshold.')
            return []

        # Create interaction matrix
        matrix = df.pivot_table(
            index=firstColumn,
            columns=secondColumn,
            values=scoreColumn,
            aggfunc='mean'
        )

        # Plot heatmap
        figWidth = max(8, len(matrix.columns) * 1.2)
        figHeight = max(6, len(matrix.index) * 0.8)

        fig, ax = plt.subplots(figsize=(figWidth, figHeight))

        heatmap = ax.imshow(
            matrix.values,
            cmap='RdYlGn',
            vmin=0,
            vmax=1,
            aspect='auto'
        )

        # Add probability values inside cells
        for i in range(len(matrix.index)):
            for j in range(len(matrix.columns)):
                value = matrix.iloc[i, j]

                if pd.notna(value):
                    ax.text(
                        j, i,
                        f'{value:.3f}',
                        ha='center',
                        va='center',
                        color='black',
                        fontsize=9
                    )

        # Axis labels
        ax.set_xticks(np.arange(len(matrix.columns)))
        ax.set_yticks(np.arange(len(matrix.index)))

        ax.set_xticklabels(
            matrix.columns,
            rotation=45,
            ha='right'
        )

        ax.set_yticklabels(matrix.index)

        ax.set_xlabel('Second structure / sequence')
        ax.set_ylabel('First structure / sequence')

        ax.set_title(
            f'SpatialPPI interaction probability heatmap\n'
            f'{filterLabel} '
            f'({len(df)} combinations)'
        )

        # Colour bar
        colorbar = fig.colorbar(heatmap, ax=ax)
        colorbar.set_label('Interaction probability')

        plt.tight_layout()
        plt.show()

        return []