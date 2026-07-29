# C++ API (preview)

Generated from the C++ sources with Doxygen and rendered by Breathe.

:::{note}
Coverage is currently **structural**. The sources carry few Doxygen
comments so far; as comments are added to the headers destined for
`libcneta`, descriptions will appear here automatically.

Regenerate locally with `make -C docs api`.
:::

## Trees

```{doxygenfile} evo_tree.hpp
:project: cneta
```

## Copy-number parsing

```{doxygenfile} parse_cn.hpp
:project: cneta
```

<!--
Add further headers as they stabilise, e.g.:

```{doxygenfile} likelihood.hpp
:project: cneta
```

Once symbols are properly documented, prefer targeted directives over
whole files:

```{doxygenclass} evo_tree
:members:
```
-->
