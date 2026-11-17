# **************************************************************************
# *
# * Authors:   Blanca Pueche (blanca.pueche@cnb.csis.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import csv
import shutil
import json
import re
import shlex

import pyworkflow.protocol.params as params

from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import cifToPdb

from pwchem import Plugin
from spatialppi import SPATIALPPI_DIC

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from pwem.objects import SetOfSequences

class ProtSpatialPPI(EMProtocol):
    """
    Predict protein-protein interactions using SpatialPPIv2.

    The protocol supports:
        - Two input structures, with optional chain selection.
        - Two individual input sequences.
        - A set of sequences for pairwise interaction screening.

    When no chains or sequences are explicitly selected in the
    corresponding screening mode, all possible combinations are analysed.
    """

    _label = 'PPI prediction'

    CSV_FILE = 'spatialppi_results.csv'
    RESULTS_FILE = 'spatialppi_results.json'
    FILTERED_CSV_FILE = 'spatialppi_filtered_results.csv'

    def _defineParams(self, form):
        form.addHidden(
            'useGpu',
            params.BooleanParam,
            default=True,
            label="Use GPU for execution",
            help="Choose whether to use the GPU or CPU."
        )

        form.addHidden(
            'gpuList',
            params.StringParam,
            default='0',
            label="Choose GPU IDs",
            help="Comma-separated GPU devices that can be used."
        )

        form.addSection(label='Input')

        form.addParam(
            'inputType',
            params.EnumParam,
            choices=['AtomStruct', 'Sequence', 'SetOfSequences'],
            label="Input type:",
            default=0,
            help='Select the input object type.'
        )

        form.addParam(
            'inputAtomStruct1',
            params.PointerParam,
            pointerClass='AtomStruct',
            condition='inputType==0',
            label='Structure 1:',
            help='Select the first structure to study.'
        )

        form.addParam(
            'chain1',
            params.StringParam,
            condition='inputType==0',
            label='Chain 1:',
            default='',
            help="Select the chains to analyse in structure 1. "
                 "Leave empty to analyse all protein chains."
        )

        form.addParam(
            'inputAtomStruct2',
            params.PointerParam,
            pointerClass='AtomStruct',
            condition='inputType==0',
            label='Structure 2:',
            help='Select the second structure to study.'
        )

        form.addParam(
            'chain2',
            params.StringParam,
            condition='inputType==0',
            label='Chain 2:',
            default='',
            help="Select the chains to analyse in structure 2. "
                 "Leave empty to analyse all protein chains."
        )

        form.addParam(
            'inputSeq1',
            params.PointerParam,
            pointerClass='Sequence',
            condition='inputType==1',
            label='Sequence 1:',
            help='Select the first sequence to study.'
        )

        form.addParam(
            'inputSeq2',
            params.PointerParam,
            pointerClass='Sequence',
            condition='inputType==1',
            label='Sequence 2:',
            help='Select the second sequence to study.'
        )

        form.addParam(
            'inputSetOfSequences',
            params.PointerParam,
            pointerClass='SetOfSequences',
            condition='inputType==2',
            label='Input set of sequences:',
            help='Select the sequence set to analyse.'
        )

        form.addParam(
            'inputSequenceFromSet1',
            params.StringParam,
            label='Sequence 1:',
            condition='inputType==2',
            default='',
            help='Select the first sequence from the set. '
                 'Leave empty to analyse all possible combinations.'
        )

        form.addParam(
            'inputSequenceFromSet2',
            params.StringParam,
            label='Sequence 2:',
            condition='inputType==2',
            default='',
            help='Select the second sequence from the set. '
                 'Leave empty to analyse all possible combinations.'
        )

        form.addParam(
            'threshold',
            params.FloatParam,
            default=0.5,
            label='Interaction probability threshold:',
            help='Minimum interaction probability required for a prediction '
                 'to pass the filter. Values range from 0 to 1.'
        )

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.runSpatialPPIStep)
        self._insertFunctionStep(self.createOutputStep)

    def runSpatialPPIStep(self):
        if self.inputType.get() == 0:
            self.runStructurePredictions()
        elif self.inputType.get() == 1:
            self.runSequencePrediction()
        else:
            self.runSequencePredictionFromSet()

    def runStructurePredictions(self):
        struct1 = self.inputAtomStruct1.get()
        struct2 = self.inputAtomStruct2.get()

        name1 = self.getStructureName(struct1)
        name2 = self.getStructureName(struct2)

        file1 = self.convertCifToPdb(
            struct1.getFileName(),
            'structure1'
        )
        file2 = self.convertCifToPdb(
            struct2.getFileName(),
            'structure2'
        )

        availableChains1 = self.getProteinChains(file1)
        availableChains2 = self.getProteinChains(file2)

        chains1 = self.selectChains(
            availableChains1,
            self.chain1.get(),
            'Structure 1'
        )
        chains2 = self.selectChains(
            availableChains2,
            self.chain2.get(),
            'Structure 2'
        )

        totalPredictions = len(chains1) * len(chains2)
        print(
            f'\nStarting SpatialPPI screening: '
            f'{name1} vs {name2}'
        )

        print(f'Structure 1 chains: {", ".join(chains1)}')
        print(f'Structure 2 chains: {", ".join(chains2)}')
        print(f'Total predictions: {totalPredictions}')

        results = []
        predictionIndex = 0

        for chain1 in chains1:
            for chain2 in chains2:
                predictionIndex += 1
                print(
                    f'\nPrediction {predictionIndex}/{totalPredictions}: '
                    f'{name1} (chain {chain1}) vs '
                    f'{name2} (chain {chain2})'
                )

                score = self.runStructurePrediction(
                    file1,
                    file2,
                    chain1,
                    chain2
                )
                print(
                    f'Interaction probability: {score:.4f}'
                )

                results.append({
                    'structure1': name1,
                    'chain1': chain1,
                    'structure2': name2,
                    'chain2': chain2,
                    'score': score
                })

        self.writeResults(results)
        print(
            f'\nScreening completed. '
            f'{len(results)} predictions were performed.'
        )

    def runStructurePrediction(self, file1, file2, chain1, chain2):
        args = [
            '--A', shlex.quote(file1),
            '--B', shlex.quote(file2),
            '--chain_A', shlex.quote(chain1),
            '--chain_B', shlex.quote(chain2),
            '--model', 'ProtT5',
            '--device', 'cuda' if self.useGpu.get() else 'cpu'
        ]

        return self.runSpatialPPICommand(args)

    def runSequencePrediction(self):
        seq1 = self.inputSeq1.get()
        seq2 = self.inputSeq2.get()

        file1 = self.convertSeqToFasta(
            seq1,
            'sequence1'
        )
        file2 = self.convertSeqToFasta(
            seq2,
            'sequence2'
        )

        score = self.runSequencePredictionPair(
            file1,
            file2
        )
        print(
            f'\nInteraction probability: {score:.4f}'
        )

        results = [{
            'protein1': seq1.getSeqName(),
            'protein2': seq2.getSeqName(),
            'score': score
        }]

        self.writeResults(results)

    def runSequencePredictionFromSet(self):
        sequenceSet = self.inputSetOfSequences.get()
        sequences = [
            seq.clone()
            for seq in sequenceSet
        ]

        selected1 = (
            self.inputSequenceFromSet1.get() or ''
        )
        selected2 = (
            self.inputSequenceFromSet2.get() or ''
        )

        seq1 = (
            self.findSequenceInSet(sequences, selected1)
        )
        print(f"seq1: {seq1}")
        seq2 = (
            self.findSequenceInSet(sequences, selected2)
        )
        print(f"seq2: {seq2}")
        if seq1 is not None and seq2 is not None:
            sequencePairs = [(seq1, seq2)]

        elif seq1 is not None:
            sequencePairs = [
                (seq1, seq)
                for seq in sequences
                if seq.getSeqName() != seq1.getSeqName()
            ]
        elif seq2 is not None:
            sequencePairs = [
                (seq, seq2)
                for seq in sequences
                if seq.getSeqName() != seq2.getSeqName()
            ]
        else:
            sequencePairs = [
                (sequences[i], sequences[j])
                for i in range(len(sequences))
                for j in range(i + 1, len(sequences))
            ]

        if not sequencePairs:
            raise RuntimeError(
                'No valid sequence combinations were found.'
            )

        totalPredictions = len(sequencePairs)

        print(
            f'\nStarting sequence screening: '
            f'{totalPredictions} prediction(s).'
        )

        results = []
        for index, (seq1, seq2) in enumerate(
            sequencePairs,
            start=1
        ):
            print(
                f'\nPrediction {index}/{totalPredictions}: '
                f'{seq1.getSeqName()} vs {seq2.getSeqName()}'
            )
            file1 = self.convertSeqToFasta(
                seq1,
                f'sequence1_{index}'
            )
            file2 = self.convertSeqToFasta(
                seq2,
                f'sequence2_{index}'
            )
            score = self.runSequencePredictionPair(
                file1,
                file2
            )
            print(
                f'Interaction probability: {score:.4f}'
            )
            results.append({
                'protein1': seq1.getSeqName(),
                'protein2': seq2.getSeqName(),
                'score': score
            })

        self.writeResults(results)
        print(
            f'\nSequence screening completed. '
            f'{len(results)} predictions were performed.'
        )

    def runSequencePredictionPair(self, file1, file2):
        args = [
            '--A', shlex.quote(file1),
            '--B', shlex.quote(file2),
            '--model', 'ESM-2+ac',
            '--device', 'cuda' if self.useGpu.get() else 'cpu'
        ]
        return self.runSpatialPPICommand(args)

    def runSpatialPPICommand(self, args):
        logFile = self._getLogsPath('run.stdout')
        previousMatches = []
        if os.path.exists(logFile):
            with open(logFile) as f:
                previousMatches = re.findall(
                    r'Possibility of interaction:\s*([0-9.eE+-]+)',
                    f.read()
                )

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=SPATIALPPI_DIC,
            program="python inference.py",
            cwd=os.path.abspath(
                os.path.join(
                    Plugin.getVar(SPATIALPPI_DIC['home']),
                    'SpatialPPIv2'
                )
            )
        )

        with open(logFile) as f:
            logContent = f.read()
        currentMatches = re.findall(
            r'Possibility of interaction:\s*([0-9.eE+-]+)',
            logContent
        )
        if not currentMatches:
            raise RuntimeError(
                'Could not retrieve the interaction probability from '
                'the SpatialPPI output. Check run.stdout.'
            )
        if len(currentMatches) > len(previousMatches):
            scoreString = currentMatches[len(previousMatches)]
        else:
            scoreString = currentMatches[-1]

        score = float(scoreString)

        if not 0.0 <= score <= 1.0:
            raise RuntimeError(
                f'SpatialPPI returned an invalid probability: {score}'
            )

        return score

    def createOutputStep(self):
        threshold = self.threshold.get()

        resultsFile = self._getExtraPath(self.RESULTS_FILE)

        with open(resultsFile) as f:
            results = json.load(f)

        filteredResults = [
            result for result in results
            if result['score'] >= threshold
        ]

        print(
            f'\nApplying interaction probability threshold: {threshold:.4f}'
        )
        print(
            f'Predictions passing the threshold: '
            f'{len(filteredResults)}/{len(results)}'
        )
        self.writeFilteredResults(filteredResults)

        if self.inputType.get() == 0:
            return

        passingSequences = set()

        for result in filteredResults:
            passingSequences.add(result['protein1'])
            passingSequences.add(result['protein2'])

        outputSequences = SetOfSequences.create(
            outputPath=self._getPath()
        )

        if self.inputType.get() == 1:
            inputSequences = [
                self.inputSeq1.get(),
                self.inputSeq2.get()
            ]
        else:
            inputSequences = self.inputSetOfSequences.get()

        addedSequences = set()

        for seq in inputSequences:
            seqName = seq.getSeqName()

            if (
                    seqName in passingSequences
                    and seqName not in addedSequences
            ):
                outputSequences.append(seq.clone())
                addedSequences.add(seqName)

        self._defineOutputs(
            outputSequences=outputSequences
        )

        if self.inputType.get() == 1:
            self._defineSourceRelation(
                self.inputSeq1,
                outputSequences
            )
            self._defineSourceRelation(
                self.inputSeq2,
                outputSequences
            )
        else:
            self._defineSourceRelation(
                self.inputSetOfSequences,
                outputSequences
            )

        print(
            f'Output sequence set created: '
            f'{outputSequences.getSize()} sequence(s).'
        )

    def writeResults(self, results):
        resultsFile = self._getExtraPath(
            self.RESULTS_FILE
        )
        with open(resultsFile, 'w') as f:
            json.dump(
                results,
                f,
                indent=4
            )

        csvFile = self._getExtraPath(
            self.CSV_FILE
        )

        if self.inputType.get() == 0:

            fieldnames = [
                'Structure 1',
                'Chain 1',
                'Structure 2',
                'Chain 2',
                'Interaction probability'
            ]

        else:

            fieldnames = [
                'Protein 1',
                'Protein 2',
                'Interaction probability'
            ]

        with open(csvFile, 'w', newline='') as f:

            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames
            )
            writer.writeheader()
            for result in results:
                if self.inputType.get() == 0:

                    writer.writerow({
                        'Structure 1': result['structure1'],
                        'Chain 1': result['chain1'],
                        'Structure 2': result['structure2'],
                        'Chain 2': result['chain2'],
                        'Interaction probability': result['score']
                    })

                else:

                    writer.writerow({
                        'Protein 1': result['protein1'],
                        'Protein 2': result['protein2'],
                        'Interaction probability': result['score']
                    })

    def writeFilteredResults(self, results):
        csvFile = self._getExtraPath(self.FILTERED_CSV_FILE)

        if self.inputType.get() == 0:
            fieldnames = [
                'Structure 1',
                'Chain 1',
                'Structure 2',
                'Chain 2',
                'Interaction probability'
            ]
        else:
            fieldnames = [
                'Protein 1',
                'Protein 2',
                'Interaction probability'
            ]

        with open(csvFile, 'w', newline='') as f:
            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames
            )

            writer.writeheader()

            for result in results:
                if self.inputType.get() == 0:
                    writer.writerow({
                        'Structure 1': result['structure1'],
                        'Chain 1': result['chain1'],
                        'Structure 2': result['structure2'],
                        'Chain 2': result['chain2'],
                        'Interaction probability': result['score']
                    })
                else:
                    writer.writerow({
                        'Protein 1': result['protein1'],
                        'Protein 2': result['protein2'],
                        'Interaction probability': result['score']
                    })

    def _summary(self):
        resultsFile = self._getExtraPath(
            self.RESULTS_FILE
        )

        if not os.path.exists(resultsFile):
            return [
                'Prediction results are not available yet.'
            ]
        with open(resultsFile) as f:
            results = json.load(f)

        if not results:
            return [
                'No predictions were performed.'
            ]

        summary = []
        if self.inputType.get() == 0:
            scores = [
                result['score']
                for result in results
            ]

            bestResult = max(
                results,
                key=lambda result: result['score']
            )

            summary.append(
                f'Structure 1: '
                f'{self.getStructureName(self.inputAtomStruct1.get())}'
            )

            summary.append(
                f'Structure 2: '
                f'{self.getStructureName(self.inputAtomStruct2.get())}'
            )

            summary.append(
                f'Number of chain combinations analysed: *{len(results)}*'
            )

            summary.append(
                f'Highest interaction probability: '
                f'*{max(scores):.4f}*'
            )

            summary.append(
                f'Best-scoring combination: '
                f'*{bestResult["chain1"]} - {bestResult["chain2"]}*'
            )

            summary.append(
                f'Lowest interaction probability: '
                f'*{min(scores):.4f}*'
            )

            summary.append(
                f'Results table: *{self.CSV_FILE}*'
            )

        elif self.inputType.get() == 1:
            result = results[0]

            summary.append(
                f'Protein 1: {result["protein1"]}'
            )

            summary.append(
                f'Protein 2: {result["protein2"]}'
            )

            summary.append(
                f'Possibility of interaction: '
                f'*{result["score"]:.4f}*'
            )

            summary.append(
                f'Results table: *{self.CSV_FILE}*'
            )
        else:
            scores = [
                result['score']
                for result in results
            ]

            bestResult = max(
                results,
                key=lambda result: result['score']
            )

            summary.append(
                f'Number of sequence combinations analysed: *{len(results)}*'
            )

            summary.append(
                f'Highest interaction probability: '
                f'*{max(scores):.4f}*'
            )

            summary.append(
                f'Best-scoring combination: '
                f'*{bestResult["protein1"]} - '
                f'{bestResult["protein2"]}*'
            )

            summary.append(
                f'Lowest interaction probability: '
                f'*{min(scores):.4f}*'
            )

            summary.append(
                f'Results table: *{self.CSV_FILE}*'
            )

        return summary

    def _validate(self):
        errors = []

        threshold = self.threshold.get()

        if threshold < 0.0 or threshold > 1.0:
            errors.append(
                'The interaction probability threshold must be '
                'between 0 and 1.'
            )

        return errors

    # ---------------------------------------------------------------------

    def _methods(self):
        return []

    def _warnings(self):
        return []

    def getStructureName(self, atomStruct):

        return os.path.splitext(
            os.path.basename(atomStruct.getFileName())
        )[0]

    # ---------------------------------------------------------------------

    def convertCifToPdb(self, inpFile, fileTag):

        pdbFile = self._getExtraPath(
            f'{fileTag}.pdb'
        )

        if inpFile.lower().endswith('.pdb'):

            if os.path.abspath(inpFile) != os.path.abspath(pdbFile):
                shutil.copy2(inpFile, pdbFile)

        else:
            if os.path.exists(pdbFile):
                os.remove(pdbFile)

            cifToPdb(
                inpFile,
                pdbFile
            )

        return os.path.abspath(pdbFile)

    # ---------------------------------------------------------------------

    def getProteinChains(self, pdbFile):

        parser = PDBParser(QUIET=True)

        structure = parser.get_structure(
            'protein',
            pdbFile
        )

        model = next(
            structure.get_models(),
            None
        )

        if model is None:
            raise RuntimeError(
                f'No models were found in {pdbFile}.'
            )

        chains = []

        for chain in model:

            hasProteinResidues = any(
                is_aa(residue, standard=False)
                for residue in chain.get_residues()
            )

            if hasProteinResidues:

                chainID = chain.id.strip()

                if not chainID:
                    raise RuntimeError(
                        f'A protein chain with a blank chain ID was found '
                        f'in {pdbFile}. SpatialPPI requires an explicit ID.'
                    )

                chains.append(chainID)

        if not chains:
            raise RuntimeError(
                f'No protein chains were found in {pdbFile}.'
            )

        return chains

    # ---------------------------------------------------------------------

    def selectChains(
        self,
        availableChains,
        selectedChains,
        structureName
    ):

        selectedChains = (
            selectedChains or ''
        ).strip()

        # Empty selection means all protein chains.
        if not selectedChains:

            print(
                f'{structureName}: no chain selected; '
                f'using all protein chains.'
            )

            return availableChains

        requestedChains = []
        try:
            chainInfo = json.loads(selectedChains)

        except (ValueError, TypeError):
            chainInfo = None

        if isinstance(chainInfo, dict):

            if 'model-chain' in chainInfo:

                # Example: {"model-chain": "0-A, 0-B"}
                for item in chainInfo['model-chain'].split(','):

                    item = item.strip()

                    if '-' not in item:
                        raise RuntimeError(
                            f'Invalid chain selection: {item}'
                        )

                    _, chainID = item.rsplit('-', 1)

                    requestedChains.append(
                        chainID.strip()
                    )

            elif 'chain' in chainInfo:
                # Example: {"model": 0, "chain": "A", ...}
                requestedChains.append(
                    str(chainInfo['chain'])
                )

            else:

                raise RuntimeError(
                    f'Unrecognised chain selection: {selectedChains}'
                )

        else:
            requestedChains = [
                chain.strip()
                for chain in selectedChains.split(',')
                if chain.strip()
            ]
        requestedChains = list(
            dict.fromkeys(requestedChains)
        )

        invalidChains = [
            chain
            for chain in requestedChains
            if chain not in availableChains
        ]

        if invalidChains:

            raise RuntimeError(
                f'Chain(s) {", ".join(invalidChains)} were not found in '
                f'{structureName}. Available protein chains: '
                f'{", ".join(availableChains)}.'
            )

        if not requestedChains:

            raise RuntimeError(
                f'No valid chains were selected for {structureName}.'
            )

        return requestedChains

    def convertSeqToFasta(self, seq, fileName):
        fastaFile = self._getExtraPath(
            f'{fileName}.fasta'
        )
        with open(fastaFile, 'w') as f:
            f.write(f'>{seq.getSeqName()}\n')
            f.write(f'{seq.getSequence()}\n')
        return os.path.abspath(fastaFile)

    # ---------------------------------------------------------------------

    def findSequenceInSet(self, sequences, selectedSequence):
        for seq in sequences:
            if seq.__str__() == selectedSequence:
                print(f'FOUND: {selectedSequence}')
                return seq.clone()
        return None