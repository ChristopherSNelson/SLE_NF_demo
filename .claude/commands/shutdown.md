# /shutdown

End-of-session procedure. Run this before closing Claude Code.

## Steps

1. Write `HANDOFF.md` in the project root with:
   - What was done this session (bullet list)
   - What's left / next steps (in priority order)
   - Key decisions made and why
   - Current blockers
   - Any active pipeline runs (check `pgrep -a -f nextflow` and note PIDs, sample sheet, outdir, expected duration)

2. Propose any additions to the Mistakes Log in `CLAUDE.md` from lessons learned this session.

3. Update `memory/MEMORY.md` with anything that should persist across sessions.

4. Stage and commit all modified files (excluding ignored paths). Use the standard co-author commit format.

5. Push to origin.

6. Print a one-line summary of the session for the user.
