import ast
import importlib
import inspect
import io
import tokenize

from docutils import nodes
from docutils.statemachine import StringList
from sphinx.util.docutils import SphinxDirective


class AutoDocGroups(SphinxDirective):
    required_arguments = 1
    has_content = False

    def run(self):
        module_name = self.arguments[0]
        module = importlib.import_module(module_name)
        filename = inspect.getsourcefile(module)

        if filename is None:
            raise self.error(f"No Python source found for {module_name}")

        self.env.note_dependency(filename)

        with tokenize.open(filename) as stream:
            source = stream.read()

        # Collect unindented DOCGROUP comments.
        events = []
        for token in tokenize.generate_tokens(io.StringIO(source).readline):
            if (
                token.type == tokenize.COMMENT
                and token.start[1] == 0
                and token.string.startswith("# DOCGROUP:")
            ):
                title = token.string.partition(":")[2].strip()
                if title:
                    events.append((token.start[0], "group", title))

        # Collect public functions and classes defined at module level.
        for item in ast.parse(source, filename=filename).body:
            if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                if "grad" in item.name or item.name.lstrip("_").startswith("derivative"):
                    continue
                #if item.name.startswith("_"): # uncomment to exclude private-members
                #    continue
       
                kind = "autoclass" if isinstance(item, ast.ClassDef) else "autofunction"
                events.append((item.lineno, kind, item.name))

        lines = []
        for _, kind, name in sorted(events):
            if kind == "group":
                lines.extend([
                    f".. rubric:: {name}",
                    "   :class: docgroup-heading",
                    "",
                ])
            else:
                lines.append(f".. {kind}:: {module_name}.{name}")
                if kind == "autoclass":
                    lines.extend([
                        "   :members:",
                        "   :member-order: bysource",
                    ])
                lines.append("")

        container = nodes.container()
        self.state.nested_parse(
            StringList(lines, source=filename),
            0,
            container,
        )
        return container.children


def setup(app):
    app.setup_extension("sphinx.ext.autodoc")
    app.add_directive("autodocgroups", AutoDocGroups)
    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
