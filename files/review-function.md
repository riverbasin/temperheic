# Review Function Against Source Paper

Review the specified function for correctness against the source literature.

For the function specified in $ARGUMENTS:

1. Read the function source code
2. Read the relevant reference card in `docs/refs/` to identify which equations apply
3. If needed, read the full PDF in `refs/` for exact equation details
4. Verify:
   - Does the implementation match the paper's equations?
   - Are sign conventions consistent (positive z downward, positive v_t downwelling)?
   - Are units consistent (SI base: meters, seconds)?
   - Is the phase convention correct (radians, atan2 quadrant handling)?
5. Check the existing tests for this function — do they cover the key cases?
6. Report: correct / issues found / suggestions

Cross-reference against `observedAmpPhase.tex` for OLS fitting equations and `Chapter_1.tex` for the broader theoretical context when relevant.
