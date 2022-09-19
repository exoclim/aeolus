# -*- coding: utf-8 -*-
"""Bells and whistles."""
from iris.experimental.representation import CubeListRepresentation


class ReprAtmoSimBase:
    """
    Produce representations of AtmoSimBase-like instance.

    This includes:
    * ``str_repr``: provides __repr__ or __str__ view as a string
    * ``html_repr``: as an HTML object, available in Jupyter notebooks.
        Specifically, this is presented as an HTML table.
    """

    _template = """
<style>
  table.aeolus {{
      white-space: pre;
      border: 1px solid;
      border-color: #f9f9ef;
      font-family: monaco, monospace;
  }}
  th.aeolus {{
      background: #084469;
      color: #fefefe;
      border-left: 1px solid;
      border-color: #0a507a;
      font-size: 1.05em;
      min-width: 50px;
      max-width: 125px;
  }}
  .aeolus-la {{
      text-align: left !important;
      white-space: pre;
  }}
  .aeolus-ra {{
      text-align: right !important;
      white-space: pre;
  }}
</style>
<table class="aeolus" id="{id}">
    {header}
    {content}
</table>
        """

    def __init__(self, atmosim):
        """
        Initialise ReprAtmoSimBase.

        Parameters
        ----------
        atmosim: aeolus.core.AtmoSimBase
            AtmoSimBase-like instance
        """
        self._cls_m = atmosim.__module__
        self._cls_n = atmosim.__class__.__name__
        self.as_id = id(atmosim)
        self.atmosim = atmosim
        self._cubes = self.atmosim._cubes
        self.name = f"{self._cls_m} {self._cls_n} '{self.atmosim.name}'"
        self._copy_attrs = [
            "description",
            "planet",
            "model",
            "model_type",
            "timestep",
            "domain",
            "vert_coord",
        ]
        self._max_len = max([len(i) for i in self._copy_attrs])
        for attr in self._copy_attrs:
            setattr(self, attr, getattr(atmosim, attr))

    def str_repr(self, short=False):
        """Represent AtmoSimBase as string."""
        summary = []
        summary.append(f"{self.name} [{len(self._cubes)} cubes]")
        # Short option - for __repr__
        if short:
            return " ".join(summary)
        # Long option - for __str__
        for attr in self._copy_attrs:
            if getattr(self, attr) is not None:
                summary.append(f"{attr:>{self._max_len}} | {getattr(self, attr)}")
        if len(self._cubes) > 0:
            summary.append("{!s}\n".format(self._cubes))
        return "\n".join(summary)

    def _make_header(self):
        top_cell = f'<th class="aeolus aeolus-la" colspan="2">{self.name}</th>'
        cells = ['<tr class="aeolus">', top_cell, "</tr>"]
        return "\n".join(cells)

    def _make_content(self):
        cells = []
        for attr in self._copy_attrs:
            if getattr(self, attr) is not None:
                value = getattr(self, attr)
                cells.append('<tr class="aeolus">')
                cells.append(f'<td class="aeolus aeolus-ra">{attr}</td>')
                cells.append(f'<td class="aeolus aeolus-la">{value!s}</td>')
                cells.append("</tr>")
        if len(self._cubes) > 0:
            # List cubes
            cl_repr_html = CubeListRepresentation(self._cubes).repr_html()
            cells.append('<tr class="aeolus">')
            cells.append('<td class="octant aeolus-ra">Cubes</td>')
            cells.append(f'<td class="aeolus aeolus-la">{cl_repr_html}</td>')
            cells.append("</tr>")

        return "\n".join(cells)

    def html_repr(self):
        """HTML representation used in Jupyter Notebooks."""
        header = self._make_header()
        content = self._make_content()
        return self._template.format(id=self.as_id, header=header, content=content)
