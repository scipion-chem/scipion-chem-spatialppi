# **************************************************************************
# *
# * Authors:  Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from .constants import *
from pwchem.constants import RDKIT_DIC

_references = ['']


class Plugin(pwchemPlugin):
    @classmethod
    def defineBinaries(cls, env):
        cls.addSpatialPPIPackage(env)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(SPATIALPPI_DIC['home'], cls.getEnvName(SPATIALPPI_DIC))

    @classmethod
    def addSpatialPPIPackage(cls, env, default=True):
        installer = InstallHelper(
            SPATIALPPI_DIC['name'],
            packageHome=cls.getVar(SPATIALPPI_DIC['home']),
            packageVersion=SPATIALPPI_DIC['version']
        )

        installer.getCondaEnvCommand(
            SPATIALPPI_DIC['name'],
            binaryVersion=SPATIALPPI_DIC['version'],
            pythonVersion='3.11'
        )

        # Basic dependencies
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "pip install numpy pandas seaborn tensorboard",
            "BASIC_DEPENDENCIES_INSTALLED"
        )

        # PyTorch Geometric
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "pip install torch_geometric==2.5.1 && "
            "pip install pyg_lib torch_scatter torch_sparse "
            "torch_cluster torch_spline_conv "
            "-f https://data.pyg.org/whl/torch-2.2.0+cu118.html",
            "PYG_INSTALLED"
        )

        # Lightning
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "conda install -y lightning -c conda-forge && "
            "pip install 'lightning[extra]'",
            "LIGHTNING_INSTALLED"
        )

        # Scientific dependencies
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "pip install "
            "matplotlib==3.8.3 "
            "biopandas==0.5.1 "
            "biopython==1.83 "
            "transformers==4.40.2 "
            "sentencepiece==0.2.0 "
            "torchsummary==1.5.1 "
            "scipy==1.12.0",
            "SCIENTIFIC_DEPENDENCIES_INSTALLED"
        )

        # Additional torch_cluster installation (as in official instructions)
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "pip install torch_cluster "
            "-f https://data.pyg.org/whl/torch-2.2.0+cu118.html",
            "TORCH_CLUSTER_INSTALLED"
        )

        # Prodigy-cryst
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "pip install git+https://github.com/yusuf1759/prodigy-cryst.git",
            "PRODIGY_INSTALLED"
        )

        # PINDER and ESM
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "pip install 'pinder[all]' && "
            "pip install fair-esm",
            "PINDER_ESM_INSTALLED"
        )

        # Clone SpatialPPIv2
        installer.addCommand(
            "git clone --depth 1 "
            "https://github.com/ohuelab/SpatialPPIv2.git",
            "SPATIALPPI_CLONED"
        )

        # Install PyTorch and CUDA LAST
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "conda install -y "
            "pytorch=2.2.0 "
            "torchvision=0.17.0 "
            "torchaudio=2.2.0 "
            "pytorch-cuda=11.8 "
            "-c pytorch -c nvidia",
            "PYTORCH_INSTALLED"
        )

        # Test imports
        installer.addCommand(
            f"{cls.getEnvActivationCommand(SPATIALPPI_DIC)} && "
            "python -c \""
            "import torch; "
            "import torchvision; "
            "import torchaudio; "
            "import torch_geometric; "
            "import torchmetrics; "
            "print('PyTorch:', torch.__version__); "
            "print('Torchvision:', torchvision.__version__); "
            "print('Torchaudio:', torchaudio.__version__); "
            "print('PyTorch Geometric:', torch_geometric.__version__); "
            "print('CUDA:', torch.version.cuda); "
            "print('CUDA available:', torch.cuda.is_available()); "
            "print('All imports OK')\"",
            "SPATIALPPI_IMPORT_TEST"
        )

        installer.addPackage(
            env,
            dependencies=['conda', 'pip', 'git'],
            default=default
        )


