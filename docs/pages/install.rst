.. _install:

.. contents::
   :depth: 3

Install
=======

To install NucKit, simply clone the repository to the desired destination.

.. code-block:: shell
  
  git clone https://github.com/cnobles/nuckit.git

Then initiate the install using the install script. If you would like the 
installed environment to be named something other than 'nuckit', the new conda 
environment name can be provided to the ``install.sh`` script as provided below.

.. code-block:: shell

  cd path/to/nuckit
  bash install.sh

Or:

.. code-block:: shell

  cd path/to/nuckit
  bash install.sh -e {env_name}
  
Additionally, help information on how to use the ``install.sh`` can be accessed
by:

.. code-block:: shell

  bash install.sh -h
