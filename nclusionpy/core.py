import os
import subprocess
from pathlib import Path

# ---------------------------------------------------------------------
# Pre-bootstrap Julia environment reproducibly
# ---------------------------------------------------------------------
repo_root = Path(__file__).resolve().parents[1]
# depot = repo_root / ".julia_depot"
# depot.mkdir(exist_ok=True)

env = dict(os.environ)
env.pop("LD_LIBRARY_PATH", None)  # avoid system curl/openssl conflicts
env["PYTHON_JULIAPKG_PROJECT"] = str(repo_root)
env["PYTHON_JULIAPKG_OFFLINE"] = "yes"


# # Enable offline mode only if already bootstrapped
# if (depot / "packages").exists():
    

cmd = (
    "using Pkg; "
    f'Pkg.activate("{repo_root}"); '
    "Pkg.resolve(); "
    'Pkg.pin(name="OpenSSL_jll", version="3.0.15"); '
    'Pkg.pin(name="PythonCall", version="0.9.24"); '
    "Pkg.resolve(); "
    "Pkg.instantiate(); "
    "Pkg.precompile();"
)

subprocess.run(
    ["julia", f"--project={repo_root}", "-e", cmd],
    check=True,
    env=env,
)

import juliapkg
import juliacall
from juliacall import Main as jl
# ---------------------------------------------------------------------
# Python-callable setup
# ---------------------------------------------------------------------
def setup_julia_environment():
    """
    Bootstraps Julia + Nclusion reproducibly on any machine (HPC or local).
    Ensures Julia is installed, activates the repo's Project.toml,
    repairs common binary mismatches, and prints a startup confirmation.
    """

    repo_root = Path(__file__).resolve().parents[1]
    depot = repo_root / ".julia_depot"

    try:
        juliapkg.resolve()
    except Exception as e:
        print("[nclusionpy] ⚠️  Failed to resolve Julia environment via juliapkg:", e)
        print("[nclusionpy] Continuing with existing Julia installation...")

    # Activate Nclusion.jl environment
    jl.seval(f'using Pkg; Pkg.activate(raw"{repo_root}")')

    # Check and fix OpenSSL_jll binary mismatch
    jl.seval("""
    try
        using Pkg
        if "OpenSSL_jll" ∈ keys(Pkg.project().dependencies)
            try
                using OpenSSL_jll
            catch e
                @warn "OpenSSL_jll failed to load, pinning to v3.0.15" exception=e
                Pkg.pin(name="OpenSSL_jll", version="3.0.15")
                Pkg.resolve()
                Pkg.instantiate()
            end
        end
        if "PythonCall" ∈ keys(Pkg.project().dependencies)
            try
                using PythonCall
            catch e
                @warn "PythonCall failed to load, pinning to v0.9.24" exception=e
                Pkg.pin(name="PythonCall", version="0.9.24")
                Pkg.resolve()
                Pkg.instantiate()
            end
        end
    catch e
        @warn "Binary dependency check failed" exception=e
    end
    """)

    # Precompile
    jl.seval("""
    using Pkg
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.precompile()
    """)

    # Load core dependencies
    jl.seval("""
    using Logging, LoggingExtras
    using JLD2, FileIO, OrderedCollections, CSV, DataFrames
    using Nclusion

    global logger = FormatLogger() do io, args
        println(io, args._module, " | ", "[", args.level, "] ", args.message)
    end
    """)

    jl_version = jl.seval("string(VERSION)")
    print(f"[nclusionpy] ✅ Julia {jl_version} detected — Nclusion.jl environment ready.")
    print(f"[nclusionpy] Activated project: {repo_root}/Project.toml\n")


def reset_julia():
    """Reset Julia Main namespace and logger."""
    jl.seval("""
    for name in names(Main)
        if !isconst(Main, name)
            Main.eval(:($name = nothing))
        end
    end
    global logger = FormatLogger() do io, args
        println(io, args._module, " | ", "[", args.level, "] ", args.message)
    end
    """)

# import os, subprocess
# from pathlib import Path
# repo_root = Path(__file__).resolve().parents[1]
# depot = repo_root / ".julia_depot"
# depot.mkdir(exist_ok=True)
# env = dict(os.environ)
# env.pop("LD_LIBRARY_PATH", None)
# env["PYTHON_JULIAPKG_PROJECT"] = str(repo_root)
# env["PYTHON_JULIAPKG_OFFLINE"] = "yes"
# env["JULIA_DEPOT_PATH"] = str(depot)
# env["JULIA_SSL_CA_ROOTS_PATH"] = "builtin"
# env["JULIA_DEPOT_PATH"] = str(repo_root / ".julia_depot")
# env["JULIA_COPY_STACKS"] = "yes"  # (optional: avoids occasional JLL init issues)
# env["JULIA_PKG_SERVER"] = ""       # disable pkg.julialang.org (use GitHub direct)
# cmd = (
#     "using Pkg; "
#     f'Pkg.activate("{repo_root}"); '
#     "Pkg.resolve(); "
#     'Pkg.pin(name="OpenSSL_jll", version="3.0.15"); '
#     "Pkg.resolve(); "
#     "Pkg.instantiate(); "
#     "Pkg.precompile();"
# )
# # Run Julia with inline code
# subprocess.run(
#     ["julia", f"--project={repo_root}", "-e", cmd],
#     check=True,
#     env=env,
# )
# import juliapkg
# from juliacall import Main as jl

# def setup_julia_environment():
#     """
#     Bootstraps Julia + Nclusion reproducibly on any machine (HPC or local).
#     Ensures Julia is installed, activates the repo's Project.toml,
#     repairs common binary mismatches, and prints a startup confirmation.
#     """

#     repo_root = Path(__file__).resolve().parents[1]

#     # === STEP 1: Ensure Julia is installed and resolve environment ===
#     try:
#         juliapkg.resolve()
#     except Exception as e:
#         print("[nclusionpy] ⚠️  Failed to resolve Julia environment via juliapkg:", e)
#         print("[nclusionpy] Continuing with existing Julia installation...")

#     # === STEP 2: Activate the repo's Project.toml ===
#     jl.seval(f"using Pkg; Pkg.activate(raw\"{repo_root}\")")

#     # === STEP 3: Fix OpenSSL_jll binary incompatibilities if needed ===
#     jl.seval("""
#     try
#         using Pkg
#         if "OpenSSL_jll" ∈ keys(Pkg.project().dependencies)
#             try
#                 using OpenSSL_jll
#             catch e
#                 @warn "OpenSSL_jll failed to load, pinning to v3.0.15" exception=e
#                 Pkg.pin(name="OpenSSL_jll", version="3.0.15")
#                 Pkg.resolve()
#                 Pkg.instantiate()
#             end
#         end
#     catch e
#         @warn "Binary dependency check failed" exception=e
#     end
#     """)

#     # === STEP 4: Resolve and precompile ===
#     jl.seval("""
#     using Pkg
#     Pkg.resolve()
#     Pkg.instantiate()
#     Pkg.precompile()
#     """)

#     # === STEP 5: Load core dependencies and logger ===
#     jl.seval("""
#     using Logging, LoggingExtras
#     using JLD2, FileIO, OrderedCollections, CSV, DataFrames
#     using Nclusion

#     global logger = FormatLogger() do io, args
#         println(io, args._module, " | ", "[", args.level, "] ", args.message)
#     end
#     """)

#     # === STEP 6: Print a confirmation banner ===
#     version_info = jl.seval("string(VERSION)")
#     jl_version = jl.seval("string(VERSION)")
#     print(f"[nclusionpy] ✅ Julia {jl_version} detected — Nclusion.jl environment ready.")
#     print(f"[nclusionpy] Activated project: {repo_root}/Project.toml\n")

# def reset_julia():
#     """Reset Main namespace and logger."""
#     jl.seval("""
#     for name in names(Main)
#         if !isconst(Main, name)
#             Main.eval(:($name = nothing))
#         end
#     end
#     global logger = FormatLogger() do io, args
#         println(io, args._module, " | ", "[", args.level, "] ", args.message)
#     end
#     """)

# import os
# # os.environ["LD_LIBRARY_PATH"] = (os.environ.get("LD_LIBRARY_PATH", "") + ":/lib64")
# from pathlib import Path
# import subprocess
# # point to your repo Project.toml
# __file__ = "/users/cnwizu/data/cnwizu/nclusion/nclusionpy/__init__.py"
# repo_root = Path(__file__).resolve().parents[1]
# os.environ["PYTHON_JULIAPKG_PROJECT"] = str(repo_root)
# os.environ["PYTHON_JULIAPKG_OFFLINE"] = "yes"

# env = dict(os.environ)
# env.pop("LD_LIBRARY_PATH", None)
# cmd = (
#     "using Pkg; "
#     f'Pkg.activate("{repo_root}"); '
#     'Pkg.pin("OpenSSL_jll", v"3.0.15+3")'
#     "Pkg.resolve(); "
#     "Pkg.instantiate(); "
#     "Pkg.precompile();"
# )
# # Run Julia with inline code
# subprocess.run(
#     ["julia", f"--project={repo_root}", "-e", cmd],
#     check=True,
#     env=env,
# )
# import juliapkg
# juliapkg.resolve(force=False, dry_run=False)

# from juliacall import Main as jl

# # Compute absolute path relative to this file
# repo_root = Path(__file__).resolve().parents[1]  # go up from nclusionpy/
# nclusion_jl = repo_root / "src" / "Nclusion.jl"


# def setup_julia_environment():
#     """
#     Ensure Julia and Nclusion dependencies are installed and available.
#     """

#     repo_root = Path(__file__).resolve().parents[1]
#     # This will install Julia itself (if missing) and instantiate Project.toml
#     juliapkg.resolve()
#     jl.seval(f'using Pkg; Pkg.activate("{repo_root}")')
#     nclusion_jl = repo_root / "src" / "Nclusion.jl"
#     jl.seval("using Pkg; Pkg.resolve()")
#     jl.seval("Pkg.instantiate()")
#     jl.seval("Pkg.precompile()")
#     jl.seval("using Nclusion")
#     jl.seval(f'include("{nclusion_jl}")')
#     # Now safe to import packages
#     jl.seval("""
#     using Logging, LoggingExtras, JLD2, FileIO, OrderedCollections, CSV, DataFrames
#     global logger = FormatLogger() do io, args
#         println(io, args._module, " | ", "[", args.level, "] ", args.message)
#     end
#     """)

# # Set up Julia environment and logger just once
# jl.seval('ENV["GKSwstype"] = "100"')
# jl.seval("""
# using Logging, LoggingExtras
# using JLD2, FileIO
# using OrderedCollections
# using CSV, DataFrames

# logger = FormatLogger() do io, args
#     println(io, args._module, " | ", "[", args.level, "] ", args.message)
# end
# """)


# # Import Julia package (already installed in their Julia environment)
# jl.seval("using Nclusion")


# # Utility: reset variables in Julia Main
# jl.seval("""
# function unbindvariables()
#     for name in names(Main)
#         if !isconst(Main, name)
#             Main.eval(:($name = nothing))
#         end
#     end
# end
# """)


# def setup_julia_environment():
#     """
#     Ensure that Nclusion and its dependencies are installed in Julia.
#     If missing, install them using Julia's Pkg API.
#     """
#     juliapkg.require_julia("1.11.1", target=None)
#     jl.seval("using Pkg")
#     # Ensure Julia environment is resolved
#     juliapkg.resolve(project=str(repo_root / "Project.toml"))
#     # Check if Nclusion is installed
#     is_installed = jl.seval("""
#         any(p -> p.name == "Nclusion", values(Pkg.project().dependencies))
#     """)

#     if not is_installed:
#         print("[nclusionpy] Installing Julia dependencies from Project.toml...")

#         # Path to your Julia Project.toml
#         repo_root = Path(__file__).resolve().parents[1]
#         project_toml = repo_root / "Project.toml"

#         if project_toml.exists():
#             jl.Pkg.activate(str(repo_root))   # activate your project
#             jl.Pkg.instantiate()              # install dependencies
#         else:
#             # fallback: just add Nclusion from registry or local dev
#             jl.Pkg.add("Nclusion")

#     # Finally, load your Julia package
#     jl.seval("using Nclusion")

#     # Set up logging etc.
#     jl.seval("""
#     using Logging, LoggingExtras, JLD2, FileIO, OrderedCollections, CSV, DataFrames

#     global logger = FormatLogger() do io, args
#         println(io, args._module, " | ", "[", args.level, "] ", args.message)
#     end
#     """)

