function serialise(svg) {
  /* Saves the figure to SVG in its current state.
   * Clones the provided SVG and sets the width/height of the clone to the
   * bounding box of the original SVG. Thus, downloaded figures will be sized
   * correctly.
   * This function returns a new Blob, which can then be downloaded.
  */
  let node = svg.node();
  const xmlns = "http://www.w3.org/2000/xmlns/";
  const xlinkns = "http://www.w3.org/1999/xlink";
  const xhtml = "http://www.w3.org/1999/xhtml";
  const svgns = "http://www.w3.org/2000/node";
  const bbox = svg.select("g").node().getBBox()

  node = node.cloneNode(true);
  node.setAttribute("width", bbox.width);
  node.setAttribute("height", bbox.height);
  node.setAttributeNS(xmlns, "xmlns", svgns);
  node.setAttributeNS(xmlns, "xmlns:xlink", xlinkns);
  node.setAttributeNS(xmlns, "xmlns:xhtml", xhtml);

  // Adjust x/y of <g> to account for axis/title position.
  // Replaces the transform attribute, so drag/zoom is ignored.
  d3.select(node)
    .select("g")
    .attr("transform", `translate(${Math.abs(bbox.x)}, ${Math.abs(bbox.y)})`)

  const serializer = new window.XMLSerializer;
  const string = serializer.serializeToString(node);
  return new Blob([string], {type: "image/node+xml"});
}

function download(blob, filename) {
  /* Downloads a given blob to filename.
   * This function appends a new anchor to the document, which points to the
   * supplied blob. The anchor.click() method is called to trigger the download,
   * then the anchor is removed.
  */
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}

function build(data) {
  const chart = new SynthasePlot.SynthasePlot()
  const plot = d3.select("#plot")
    .selectAll("div")
    .data([data])
    .join("div")
    .attr("class", "synthaserPlot")

  plot.call(chart)

  // General plot settings
  d3.select("#input-plotwidth")
    .on("input", function() { chart.plotWidth(+this.value); plot.call(chart) })
  d3.select("#input-title-fontsize")
    .on("input", function() { chart.titleFontSize(+this.value); plot.call(chart) })
  d3.select("#input-xlabel-fontsize")
    .on("input", function() { chart.xLabelFontSize(+this.value); plot.call(chart) })

  // Synthase bars
  d3.select("#input-barheight")
    .on("input", function() { chart.barHeight(+this.value); plot.call(chart) })
  d3.select("#input-headwidth")
    .on("input", function() { chart.headWidth(+this.value); plot.call(chart) })
  d3.select("#input-ylabel-fontsize")
    .on("input", function() { chart.yLabelFontSize(+this.value); plot.call(chart) })
  d3.select("#input-ylabel-gap")
    .on("input", function() { chart.yAxisGap(+this.value); plot.call(chart) })
  d3.select("#input-ylabel-padding")
    .on("input", function() { chart.yAxisPadding(+this.value); plot.call(chart) })

  // Legend
  d3.select("#input-legend-cellheight")
    .on("input", function() { chart.cellHeight(+this.value); plot.call(chart) })
  d3.select("#input-legend-cellwidth")
    .on("input", function() { chart.cellWidth(+this.value); plot.call(chart) })
  d3.select("#input-legend-cellpadding")
    .on("input", function() { chart.cellPadding(+this.value); plot.call(chart) })
  d3.select("#input-legend-fontsize")
    .on("input", function() { chart.legendFontSize(+this.value); plot.call(chart) })

  // Save SVG file by clicking button
  let svg = plot.select("svg")
  d3.select("#btn-save-svg")
    .on("click", () => {
      const blob = serialise(svg)
      download(blob, "synthaser.svg")
    })
}

if (typeof data === 'undefined') {
  d3.json("data.json").then(build)
} else {
  build(data)
}
