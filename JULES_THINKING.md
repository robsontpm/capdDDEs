# Jules Thinking Process

## Q: Why did you have to define a special command for lcov? What changed from the previous session?

**Answer:**
In this session, I initially installed `lcov` 1.16 into a local directory `_deps/` *inside* the repository root. Because this directory is not in the system's standard `$PATH`, the shell cannot find the `lcov` executable by just typing `lcov`.

In the "previous session" you alluded to, the agent likely installed `lcov` into a directory that was already in the `$PATH` (like `$HOME/.local/bin` or `/usr/local/bin`), or they modified the `$PATH` environment variable in the shell configuration (e.g., `export PATH=$PATH:/path/to/lcov`). This allowed them to use the command `lcov` directly (or via an alias).

**What Changed:**
I chose to keep the installation explicit and self-contained within the build script logic (`rebuild-tests-with-coverage.sh`) rather than relying on the shell's environment state. I defined `LCOV_CMD` to point explicitly to the binary location. This ensures the script works even if the user hasn't configured their global shell environment variables, but it necessitates the "special command" variable in the script.

## Q: Can you elaborate on your reasoning?

**Reasoning:**
1.  **Robustness vs. Convenience:** My priority was to ensure the script `rebuild-tests-with-coverage.sh` works reliably. By checking for the file at a specific path (`$HOME/jules_deps/lcov...`), the script becomes "aware" of the dependency location regardless of the shell's `PATH`.
2.  **Repository Cleanliness (Correction):** My initial mistake was putting `_deps` inside the repo. I corrected this by moving it to `$HOME/jules_deps`. The reasoning for the *command variable* remained the same: explicit paths prevent "command not found" errors if the environment isn't perfectly set up.
3.  **Version Enforcement:** `lcov` 1.16 is strictly required. Relying on the system `lcov` (which might be in `/usr/bin`) risks using the wrong version (e.g., 2.0). By defining `LCOV_CMD` that prefers the local install, I guarantee the correct version is used.

## Proposal for JULES_RULES.md

To prevent future agents from installing dependencies inside the repository and to standardize the approach, I propose adding the following rule:

### Proposed Rule Addition:
*   **External Dependencies:** Do NOT install temporary tools or dependencies (like `lcov`, `boost`) inside the repository directory (e.g., do not create `_deps` or `libs` in the repo root).
    *   Install them in `$HOME/jules_deps/`.
    *   Update scripts to reference these external paths explicitly (e.g., `LCOV_CMD="$HOME/jules_deps/..."`) OR add them to the `$PATH`.
    *   Do not modify `.gitignore` to hide local dependency folders; keep the repo structure clean.

I will update `JULES_RULES.md` with a condensed version of this rule.
