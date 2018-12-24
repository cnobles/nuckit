#!/usr/bin/env bash

__conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

read -r -d '' __usage <<-'EOF'
  -e --environment  [arg] Environment to install to. Default: "nuckit"
  -n --nuckit_dir   [arg] Location of NucKit source code. Default: this directory
  -c --conda  [arg]       Location of Conda installation. Default: ${PREFIX}
  -u --update [arg]       Update NucKit [tools], conda [env], or [all].
  -w --without_conda      Install without setting up conda or environment.
  -i --ignore_ctrl_mod    Ignore installation of control module.
  -r --cran               Install missing dependencies from CRAN.
  -m --cran_mirror  [arg] Install missing dependencies from supplied mirror URL.
  -v --verbose            Show subcommand output
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

read -r -d '' __helptext <<-'EOF'
 This script installs or upgrades NucKit, including Conda (if not installed).
 To upgrade, pass the '--upgrade all' option.
EOF

# Load BASH3Boilerplate for command-line parsing and logging
source etc/b3bp.sh

function __err_report() {
    local error_code
    error_code=${?}
    error "Error in ${__file} in function ${1} on line ${2}"
    exit ${error_code}
}
trap '__err_report "${FUNCNAME:-.}" ${LINENO}' ERR

# help mode
if [[ "${arg_h:?}" = "1" ]]; then
    # Help exists with code 1
    help "Help using ${0}:"
fi

# verbose mode
if [[ "${arg_v:?}" = "1" ]]; then
    LOG_LEVEL="7"
fi

# debug mode
if [[ "${arg_d:?}" = "1" ]]; then
    set -o xtrace
    LOG_LEVEL="7"
fi

function debug_capture () {
    debug "$(echo -e "$(${@})")"
}

function installation_error () {
    error "${1} failed!"
    if [[ "${arg_v:?}" != 1 && "${arg_d:?}" != 1 ]]; then
        error "Try re-running with -v or -d, or file an issue on Github."
    fi
    exit 1
}

# Set variables
__conda_path="${arg_c:-${HOME}/miniconda3}"
__nuckit_dir="${arg_n:-$(readlink -f ${__dir})}"
__nuckit_env="${arg_e:-nuckit}"
__update_tools=false
__update_env=false
__with_conda=true
__install_ctrl_mod=true
__install_cran=false
__install_mirror=false
__req_r_version="3.4.1"
__old_path=$PATH

# update env and/or tools
if [[ "${arg_u}" = "all" || "${arg_u}" = "env" ]]; then
    __update_tools=true
    __update_env=true
elif [[ "${arg_u}" = "tools" ]]; then
    __update_tools=true
fi

# install without conda
if [[ "${arg_w:?}" == "1" ]]; then
    __with_conda=false
fi

# ignore installation of control module
if [[ "${arg_i:?}" == "1" ]]; then
    __install_ctrl_mod=false
fi

# install mode for missing dependencies from cran or a cran mirror
if [[ "${arg_r:?}" == "1" || "${arg_m:-1}" == "1" ]]; then
    __install_cran=true
    __install_mirror=${arg_m:-false}
fi

# set Path
if [[ $__with_conda == true ]]; then
    PATH=$PATH:${__conda_path}/bin
fi


function __test_conda () {
    command -v conda &> /dev/null && echo true || echo false
}

function __detect_conda_install () {
    local discovered=$(__test_conda)
    if [[ $discovered = true ]]; then
        local conda_path="$(which conda)"
        echo ${conda_path%'/bin/conda'}
    fi
}

function __test_r_version () {
    local sem_version=$(R --version | grep 'R version' | cut -d ' ' -f 3)
    if (( ${sem_version//./} >= ${__req_r_version//./} )); then
        echo true
    else
        echo false
    fi
}

function __test_r_packages () {
    if [[ $__with_conda == true ]]; then activate_nuckit; fi
    $(Rscript ${__nuckit_dir}/etc/check_for_required_packages.R > /dev/null) && \
        echo true || echo false
    if [[ $__with_conda == true ]]; then deactivate_nuckit; fi
}

function __test_env () {
    if [[ $(__test_conda) = true ]]; then
          $(conda env list \
                | cut -d ' ' -f 1 \
                | grep -Fxq $__nuckit_env > /dev/null) && \
          echo true || echo false
    else
          echo false
    fi
}

function __test_nuckit () {
    if [[ $(__test_env) = true ]]; then
          activate_nuckit
          command -v nuc &> /dev/null && echo true || echo false
          deactivate_nuckit
    elif [[ $__with_conda == false ]]; then
        command -v nuc &> /dev/null && echo true || echo false
    else
          echo false
    fi
}

function activate_nuckit () {
    set +o nounset
    source activate $__nuckit_env
    set -o nounset
}

function deactivate_nuckit () {
    set +o nounset
    source deactivate
    set -o nounset
}

function install_conda () {
    local tmpdir=$(mktemp -d)
    debug "Downloading miniconda..."
    debug_capture wget -nv ${__conda_url} -O ${tmpdir}/miniconda.sh 2>&1
    debug "Installing miniconda..."
    debug_capture bash ${tmpdir}/miniconda.sh -b -p ${__conda_path} 2>&1
    if [[ $(__test_conda) != true ]]; then
        installation_error "Environment creation"
    fi
    rm ${tmpdir}/miniconda.sh
}

function install_environment () {
    debug_capture conda env update --name=$__nuckit_env \
              --quiet --file etc/environment.yml
    if [[ $(__test_env) != true ]]; then
          installation_error "Environment creation"
    fi
}

function install_env_vars () {
    if [[ $__with_conda == true ]]; then
        activate_nuckit
        echo -ne "#/bin/sh\nexport NUCKIT_DIR=${__nuckit_dir}" > \
              ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
        echo -ne "#/bin/sh\nunset NUCKIT_DIR" > \
              ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh
        deactivate_nuckit
    else
        echo -ne "\nexport NUCKIT_DIR=${__nuckit_dir}" >> ~/.bashrc
    fi
}

function install_ctrl_mod_lib () {
    if [[ $__with_conda == true ]]; then
        activate_nuckit
        debug_capture pip install --upgrade ${__nuckit_dir}/nuckit_ctrl 2>&1
        deactivate_nuckit
    else
        debug_capture pip install --upgrade ${__nuckit_dir}/nuckit_ctrl 2>&1
    fi

    if [[ $(__test_nuckit) != true ]]; then
        installation_error "Control module installation"
    fi
}

function install_nuckit_requirements () {
    if [[ $(__test_r_version) == false ]]; then
        installation_error "Insufficient R-program version"
    fi

    local check_args=" -q"

    if [[ $__install_cran == true ]]; then
        if [[ $__install_mirror != false ]]; then
            check_args=$(echo ${check_args} " --cran_mirror $__install_mirror")
        else
            check_args=$(echo ${check_args} " --cran")
        fi
    fi

    if [[ $__with_conda == true ]]; then
        check_args=$(echo ${check_args} "--conda")
    fi

    debug_capture Rscript ${__nuckit_dir}/etc/check_for_required_packages.R \
        ${check_args} 2>&1 

    if [[ $(__test_r_packages) == false ]]; then
        installation_error "R-package installation"
    fi
}

function test_nuckit_tools () {
    if [[ $__with_conda == true ]]; then activate_nuckit; fi

    if [[ $(__test_nuckit) == true ]]; then
        debug_capture bash ${__nuckit_dir}/etc/test_ctrl.sh 2>&1
    else
        debug_capture bash ${__nuckit_dir}/etc/test.sh 2>&1
    fi

    if [[ $__with_conda == true ]]; then deactivate_nuckit; fi
}


info "Starting NucKit installation..."
if [[ $__with_conda == true ]]; then
    info "    Conda path:   ${__conda_path}"
    info "    NucKit env:   ${__nuckit_env}"
fi
info "    NucKit src:   ${__nuckit_dir}"

debug "Components detected:"
if [[ $__with_conda == true ]]; then

    __conda_installed=$(__test_conda)
    debug "    Conda:               ${__conda_installed}"
    __env_exists=$(__test_env)
    debug "    Environment:         ${__env_exists}"
    
    if [[ $__env_exists == true ]]; then
        activate_nuckit
        __nuckit_installed=$(__test_nuckit)
        deactivate_nuckit
    else
        __nuckit_installed=$(__test_nuckit)
    fi
    
    debug "    Control Module:     ${__nuckit_installed}"

    __env_changed=false

    # Install Conda if necessary
    if [[ $__conda_installed = true ]]; then
        if [[ $(__detect_conda_install) != $__conda_path ]]; then
              warning "Found pre-existing Conda installation in $(__detect_conda_install)".
              warning "Ignoring specified Conda path in favor of existing Conda install."
              __conda_path=$(__detect_conda_install)
        fi
        info "Conda already installed."
    else
        info "Installing Conda..."
        install_conda
        __env_changed=true
    fi

    # Create Conda environment for NucKit
    if [[ $__env_exists = true && $__update_env = false ]]; then
        info "Specified environment already exists (use '--update env' to update)"
    else
        info "Creating NucKit environment..."
        install_environment
        __env_changed=true
    fi

    # Check environment for required r-packages
    info "Checking R-package requirements..."
    activate_nuckit
    install_nuckit_requirements
    deactivate_nuckit
    
else

    __r_installed=$(__test_r_version)
    debug "    R program:          ${__r_installed}"
    __nuckit_installed=$(__test_nuckit)
    debug "    Control Module:     ${__nuckit_installed}"
    __env_changed=false

    # Check / install R requirements
    info "Checking/Installing R-package requirments..."
    install_nuckit_requirements

fi

# Install ctrl module into environment if changed or requested
if [[ $__install_ctrl_mod == true ]]; then
    if [[ $__env_changed = true ]]; then
        info "Environment installed/updated; (re)installing NucKit library..."
        install_ctrl_mod_lib
    elif [[ $__nuckit_installed = false ]]; then
        info "Installing NucKit control module library..."
        install_ctrl_mod_lib
    elif [[ $__update_tools = true ]]; then
        info "Updating NucKit control module library..."
        install_ctrl_mod_lib
    else
        info "NucKit control module library already installed (use '--update tools' to update)"
    fi
fi

# Always update the env_vars.sh in the nuckit environment
debug "Updating \$NUCKIT_DIR variable to point to ${__nuckit_dir}"
install_env_vars

# Test NucKit tools for functionality
info "Testing NucKit tools..."
test_nuckit_tools

# Check if on pre-existing path
if [[ $__old_path != *"${__conda_path}/bin"* ]]; then
    warning "** Conda was not detected on your PATH. **"
    warning "This is normal if you haven't installed Conda before."
    warning "To add it to your path, run:"
    warning "   'echo \"export PATH=\$PATH:${__conda_path}/bin\" >> ~/.bashrc'"
    warning "and close and re-open your terminal session to apply."
    warning "When finished, run 'source activate ${__nuckit_env}' to begin."
else
    info "Done. Run 'source activate ${__nuckit_env}' to begin."
fi
