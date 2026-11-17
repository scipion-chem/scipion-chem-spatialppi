# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.tests import setupTestProject, DataSet, BaseTest

# Scipion chem imports
from pwchem.protocols import ProtChemImportSmallMolecules, ProtDefineStructROIs
from pwem.protocols import ProtImportPdb, ProtImportSequence
from pwchem.utils import assertHandle

from ..protocols import ProtSpatialPPI

class TestSpatialPPI(BaseTest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)

        cls.str1 = cls._runImportPDB('5ni1.pdb')
        cls.str2 = cls._runImportPDB('3lqd.cif')

        cls.seq1 = cls._runImportSeq('3lqd_B_mutated.fasta')
        cls.seq2 = cls._runImportSeq('AJ520101.fasta')

        cls._waitOutput(cls.str1, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.str2, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.seq1, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.seq2, 'outputSequence', sleepTime=5)

    @classmethod
    def _runImportPDB(cls, name):
        protImportPDB = cls.newProtocol(
            ProtImportPdb, inputPdbData=1,
            pdbFile=cls.ds.getFile(f'PDBx_mmCIF/{name}'))
        cls.proj.launchProtocol(protImportPDB, wait=True)
        return protImportPDB

    @classmethod
    def _runImportSeq(cls, name):
        protImportSeq = cls.newProtocol(
            ProtImportSequence, inputSequenceName=name,
            inputProteinSequence=2,
            fileSequence=cls.ds.getFile(f'Sequences/{name}'))
        cls.proj.launchProtocol(protImportSeq, wait=True)
        return protImportSeq

    def _runSpatialPPI(self, choice, obj1, obj2):
        if choice == 0:
            protPPI = self.newProtocol(ProtSpatialPPI,
                                       useGpu=False,
                                       inputType=choice,
                                       inputAtomStruct1=obj1,
                                       chain1='A',
                                       chain2='A',
                                       inputAtomStruct2=obj2
                                       )
        elif choice == 1:
            protPPI = self.newProtocol(ProtSpatialPPI,
                                       useGpu=False,
                                       inputType=choice,
                                       inputSeq1=obj1,
                                       inputSeq2=obj2
                                       )

        self.proj.launchProtocol(protPPI, wait=True)
        return protPPI

    def _waitSummary(self, protocol, expected, timeout=600, sleepTime=10):
        import time
        elapsed = 0
        while elapsed < timeout:
            summary = protocol._summary()
            if any(expected.lower() in line.lower() for line in summary):
                return
            time.sleep(sleepTime)
            elapsed += sleepTime

        self.fail(f"'{expected}' not found in protocol summary after {timeout} seconds. "
                  f"Summary: {protocol._summary()}")

    def test(self):
        protPPIstr = self._runSpatialPPI(0, self.str1.outputPdb, self.str2.outputPdb)
        self._waitSummary(protPPIstr, 'results table')

        summary = protPPIstr._summary()
        print(summary)
        self.assertTrue(any('results table' in line.lower() for line in summary))

        protPPIseq = self._runSpatialPPI(1, self.seq1.outputSequence, self.seq2.outputSequence)
        self._waitOutput(protPPIseq, 'outputSequences', sleepTime=5)
        all = getattr(protPPIseq, 'outputSequences', None)
        self.assertIsNotNone(all)


