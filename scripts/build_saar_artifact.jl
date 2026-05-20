#####
##### Build the SAAR-data artifact tarball + hashes for upload to NumericalEarth/NumericalEarthArtifacts.
#####
##### Usage:   julia --project=. scripts/build_saar_artifact.jl
#####
##### Produces `build/gsw_saar_data.tar.gz` and prints the `git-tree-sha1` and `sha256` to paste into
##### `Artifacts.toml` once the tarball is uploaded.
#####

using Pkg.Artifacts: create_artifact, archive_artifact
using Pkg.GitTools

const PKG_ROOT     = normpath(joinpath(@__DIR__, ".."))
const SOURCE_FILE  = joinpath(PKG_ROOT, "data", "gsw_saar_data.bin")
const BUILD_DIR    = joinpath(PKG_ROOT, "build")
const TARBALL_PATH = joinpath(BUILD_DIR, "gsw_saar_data.tar.gz")

isfile(SOURCE_FILE) || error("Source file not found at $SOURCE_FILE")
mkpath(BUILD_DIR)

# Stage the file inside an artifact directory; the resulting `git-tree-sha1` is the artifact's identifier.
tree_hash = create_artifact() do dir
    cp(SOURCE_FILE, joinpath(dir, "gsw_saar_data.bin"))
end

# Tar+gzip and compute the sha256 of the tarball (used by `Artifacts.toml` to verify the download).
tarball_sha256 = archive_artifact(tree_hash, TARBALL_PATH)

println()
println("Artifact built:")
println("    tarball         = ", TARBALL_PATH)
println("    git-tree-sha1   = ", tree_hash)
println("    sha256          = ", tarball_sha256)
println()
println("Paste the following stanza into `Artifacts.toml` once the tarball is uploaded, and replace")
println("`<URL>` with the public download URL on NumericalEarth/NumericalEarthArtifacts:")
println()
println("[gsw_saar_data]")
println("git-tree-sha1 = \"", tree_hash, "\"")
println("lazy = true")
println()
println("    [[gsw_saar_data.download]]")
println("    sha256 = \"", tarball_sha256, "\"")
println("    url    = \"<URL>\"")
