{pkgs ? import ./nix/pkgs.nix {}}:
pkgs.stdenv.mkDerivation {
    name = "adrenaline";
    src = pkgs.lib.cleanSource ./.;
    buildInputs = [
        pkgs.cargo
    ];
    buildPhase = ''
        export RUSTFLAGS='--emit asm'
        cargo test --release --color always
        cargo build --release --color always
        cargo doc --release --color always
    '';
    installPhase = ''
        mkdir -p "$out/bin" "$out/share"
        mv 'target/release/adrenaline' "$out/bin"
        mv 'target/doc' "$out/share"
        # TODO: Compress the assembly dumps.
        mkdir "$out/share/debug" &&                                         \
            mv 'target/release/deps/'*'.s' "$out/share/debug"
    '';
}
