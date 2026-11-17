=======================
SpatialPPIv2 plugin
=======================

**Documentation under development, sorry for the inconvenience**

This is a **Scipion** plugin that offers `SpatialPPIv2 <https://github.com/ohuelab/SpatialPPIv2.git>`_


==========================
Install this plugin
==========================

You will need to first install
`Scipion3 <https://scipion-em.github.io/docs/release-3.0.0/docs/scipion-modes/how-to-install.html>`_  and
`Scipion-chem <https://github.com/scipion-chem/scipion-chem>`_ to run these protocols.


1. **Install the plugin in Scipion**

DrugCLIP is installed automatically by scipion.

- **Install the stable version (Not available yet)**

    Through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

    or

.. code-block::

    scipion3 installp -p scipion-chem-spatialppi


- **Developer's version**

    1. **Download repository**:

    .. code-block::

        git clone https://github.com/scipion-chem/scipion-chem-spatialppi.git

    2. **Switch to the desired branch** (master or devel):

    Scipion-chem-spatialppi is constantly under development and including new features.
    If you want a relatively older an more stable version, use master branch (default).
    If you want the latest changes and developments, user devel branch.

    .. code-block::

                cd scipion-chem-spatialppi
                git checkout devel

    3. **Install**:

    .. code-block::

        scipion3 installp -p path_to_scipion-chem-spatialppi --devel




