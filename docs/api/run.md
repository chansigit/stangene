# stangene (top-level)

The top-level `stangene` module provides the main pipeline entry point and re-exports all public functions.

## Pipeline

```{eval-rst}
.. autofunction:: stangene.run
```

## Re-exported functions

All public functions are available directly from `stangene`:

```python
import stangene

stangene.run(...)
stangene.build_reference(...)
stangene.load_reference(...)
stangene.load_features(...)
stangene.classify_features(...)
stangene.harmonize(...)
stangene.merge_features(...)
stangene.summary(...)
stangene.conflict_report(...)
stangene.generate_markdown_report(...)
stangene.write_reports(...)
```
