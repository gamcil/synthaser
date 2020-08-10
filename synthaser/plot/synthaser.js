/* synthaser plot
 * Cameron L.M. Gilchrist, 2020.
 *
 * TODO: rename?
 * a utomatic
 * d omain
 * a rchitecture
 * p rofiling
 * t ool
 */

const constants = {
	"maxWidth": 600,
	"barHeight": 20,
	"headWidth": 10,
  "cellWidth": 30,
	"cellHeight": 20,
	"cellPadding": 0.3,
	"scaleFactor": 5
}

if (typeof data === 'undefined') {
	const data = d3.json("data.json").then(data => plot(data))
} else {
	plot(data);
}


function removeProps(obj, props) {
	/* Returns a copy of an object without specific properties.
	*/
	return Object.fromEntries(
		Object.keys(obj)
			.filter(p => !props.includes(p))
			.map(p => [ p, obj[p] ])
	)
}

function getUniqueDomainProps(data, prop="type") {
	/* Get all unique property values in sequence domains.
	*/
	return [
		...data.order.map(h => {
			let doms = data.synthases[h].domains.map(d => d[prop])
			return new Set(doms)
		}).reduce((one, two) => new Set([...one, ...two]))
	]
}

function getSummaryHTML(data) {
	/* Generates the HTML content for the count summary.
	 */
	let doms = getUniqueDomainProps(data, "type").length
	let fams = getUniqueDomainProps(data, "domain").length
	return `<p>Your search of <b>${data.order.length}</b> sequences`
		+ ` identified <b>${doms}</b> unique domain types`
		+ ` represented by <b>${fams}</b> different conserved domain families.`
}

function serialise(svg) {
	/* Saves the figure to SVG in its current state.
	 * Clones the provided SVG and sets the width/height of the clone to the
	 * bounding box of the original SVG. Thus, downloaded figures will be sized
	 * correctly.
	 * This function returns a new Blob, which can then be downloaded.
	*/
	node = svg.node();
	const xmlns = "http://www.w3.org/2000/xmlns/";
	const xlinkns = "http://www.w3.org/1999/xlink";
	const svgns = "http://www.w3.org/2000/node";
	const bbox = svg.select("g").node().getBBox()

	node = node.cloneNode(true);
	node.setAttribute("width", bbox.width);
	node.setAttribute("height", bbox.height);
	node.setAttributeNS(xmlns, "xmlns", svgns);
	node.setAttributeNS(xmlns, "xmlns:xlink", xlinkns);

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

function copyToClipboard(text) {
	let dummy = document.createElement("textarea")
	document.body.appendChild(dummy)
	dummy.value = text
	dummy.select()
	document.execCommand("copy")
	document.body.removeChild(dummy)
}

function getTooltipHTML(d, data) {
	/* Generates the HTML content for a cell hovering tooltip.
	 * It provides the name of the query, as well as a table of each hit.
	 */
	let pSequence = data.synthases[d.parent].sequence
	let dSequence = pSequence.slice(d.start, d.end)
	let cddUrl = "https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="
	return `
	<p class="tooltip-summary">
		<span><b>${d.parent}: ${d.type}</b></span>
	</p>
	<table class="tooltip-hits">
	<tbody>
		<tr>
			<td><b>Family</b></td>
			<td><a href="${cddUrl}${d.accession}">${d.domain}</a></td>
		</tr>
		<tr>
			<td><b>Superfamily</b></td>
			<td><a href="${cddUrl}${d.superfamily}">${d.superfamily}</a></td>
		</tr>
		<tr>
			<td><b>Class</b></td>
			<td>${d.type}</td>
		</tr>
		<tr>
			<td><b>Position</b></td>
			<td>${d.start}-${d.end}</td>
		</tr>
		<tr>
			<td><b>E-value</b></td>
			<td>${d.evalue}</td>
		</tr>
		<tr>
			<td><b>Bitscore</b></td>
			<td>${d.bitscore}</td>
		</tr>
		<tr>
		<td colspan=2>
			<b>Copy sequence:</b>
			<br>
			<button style="margin-bottom: 2px" onClick="copyToClipboard('${dSequence}')">
				Domain (${dSequence.length}aa)
			</button>
			<br>
			<button onClick="copyToClipboard('${pSequence}')">
				Protein (${pSequence.length}aa)
			</button>
		</td>
		</tr>
	</tbody>
	</table>
	`
}

function plot(data) {
	const plotDiv = d3.select("#plot");
	const svg = plotDiv.append("svg")
		.classed("wrapper-svg", true)
		.attr("id", "root_svg")
		.attr("cursor", "grab")
		.attr("xmlns", "http://www.w3.org/2000/svg")
	const transform = d3.zoomIdentity
		.translate(300, 50)
		.scale(1.2)

	// Set up the save SVG button
	d3.select("#btn-save-svg")
		.on("click", () => {
			const blob = serialise(svg)
			download(blob, "synthaser.svg")
		})

	// Reset to the original data. Have to make a deep copy here, since update
	// will mutate data
	d3.select("#btn-reset-filters")
		.on("click", () => {
			const copy = JSON.parse(JSON.stringify(data))
			update(copy)
		})

	const g = svg.append("g")
	const title = g.append("text")
		.attr("class", "title")
		.attr("text-anchor", "middle")
		.style("font-size", 20)
		.style("font-weight", "bold")
	const zoom = d3.zoom()
		.scaleExtent([0, 8])
		.on("zoom", () => g.attr("transform", d3.event.transform))
		.on("start", () => svg.attr("cursor", "grabbing"))
		.on("end", () => svg.attr("cursor", "grab"))

	svg.call(zoom).call(zoom.transform, transform)

	const x = d3.scaleLinear()
	const y = d3.scaleBand().padding(0.4)

	// Figure skeleton
	const chart = g.append("g").attr("transform", "translate(0, 10)")
	const gLegend = g.append("g")
	const bar = chart.append("g")
	const yAxis = chart.append("g")
	const xAxis = chart.append("g")
	const xLabel = chart.append("text")
		.text("Sequence length (amino acids)")
		.attr("text-anchor", "middle")
		.style("font-size", "14px")
	const annotations = chart.append("g")
	const picker = plotDiv.append("input")
		.attr("id", "picker")
		.attr("type", "color")
		.style("opacity", 0)
	const defs = svg.append("defs")

	// Tell each domain its parent
	data.order.forEach(s => {
		data.synthases[s].domains.forEach(d => {
			d["parent"] = s
			d["pLength"] = data.synthases[s].sequence.length
		})
	})

	// Set up scales for the legend based on original data domains
	const domains = getUniqueDomainProps(data)
	const scheme = d3.scaleSequential(d3.interpolateRainbow)
	const colors = d3.scaleOrdinal()
		.domain(domains)
		.range(domains.map((_, i) => scheme(i / domains.length)))
	const order = d3.scaleOrdinal()
		.domain(domains)
		.range(domains.map((_, i) => i / domains.length))

	const calculatePoints = (d) => {
		/* Function for calculating points of each sequence polygon.
		*/
		let wide = x(d.sequence.length)
		let head = wide - constants.headWidth
		let path = `
			0,0
			${head},0
			${wide},${y.bandwidth() / 2}
			${head},${y.bandwidth()}
			0,${y.bandwidth()}
		`
		return path
	}

	const cellEnter = (d, i, n) => {
		/* Populates tooltip with current cell data, and adjusts position to match 
		 * the cell in the heatmap (ignoring <g> transforms).
		*/
		tooltip.html(getTooltipHTML(d, data))
		let rect = n[i].getBoundingClientRect()
		let bbox = tooltip.node().getBoundingClientRect()
		let xOffset = rect.width / 2 - bbox.width / 2
		let yOffset = rect.height * 1.2
		tooltip
			.style("left", rect.x + xOffset + "px")
			.style("top", rect.y + yOffset + "px")
		tooltip.transition()
			.duration(100)
			.style("opacity", 1)
			.style("pointer-events", "all")
	}

	const tooltipEnter = () => {
		/* Transition upon entering the tooltip <div>.
		 * This cancels out a previously called transition (i.e. delayed transition
		 * in tooltipLeave). Also enables pointer events to allow text selection,
		 * clicking hyperlinks, etc.
		 */
		tooltip.transition()
			.duration(0)
			.style("opacity", 1)
			.style("pointer-events", "all")
	}

	const tooltipLeave = () => {
		/* Delayed tooltip transition for either 1) when user has left heatmap cell
		 * and does not go into another, or 2) user entered tooltip <div> and has now
		 * left it. Hides tooltip after 400ms, and disables pointer events which would
		 * swallow pan/zoom events.
		 */
		tooltip.transition()
			.delay(400)
			.style("opacity", 0)
			.style("pointer-events", "none")
	}

	const tooltip = plotDiv.append("div")
		.classed("tooltip", true)
		.style("opacity", 0)
		.style("pointer-events", "none")
		.style("position", "absolute")
		.style("padding", "5px")
		.on("mouseenter", tooltipEnter)
		.on("mouseleave", tooltipLeave)

	function update(data) {
		/* Updates the plot with new data.
		*/
		const t = d3.transition().duration(200)

		// Populate the summary
		d3.select("#p-result-summary")
			.html(getSummaryHTML(data))

		// Adjust range given new width/height parameters
		let maxSeq = d3.max(data.order.map(h => data.synthases[h].sequence.length))
		console.log("maxwidth", constants.maxWidth)
		x.domain([0, maxSeq])
			.range([0, maxSeq / constants.scaleFactor]) //constants.maxWidth])
			// .nice()
		y.domain(data.order)
			.range([0, data.order.length * constants.barHeight])

		const points = {}
		data.order.forEach(header => {
			points[header] = calculatePoints(data.synthases[header])
		})

		// Add <clipPath> elements for each synthase
		defs.selectAll("clipPath")
			.data(data.order, s => `${s}-clip-group`)
			.join(
				enter => enter.append("clipPath")
					.attr("id", s => `${s}-clip`)
					.append("polygon")
					.attr("points", s => points[s]),
				update => update.call(
					update => {
						return update.selectAll("polygon")
							.transition(t)
							.attr("points", s => points[s])
					}
				)
			)

		// Add <g> elements for each synthase, containing a background <rect>
		// and inner <g> for domain <rect> elements.
		const translate = (d) => `translate(0, ${y(d)})`

		// Draw the sequence bar groups.
		// Each consists of an underlying white background rect element and group of
		// domain rect elements, clipped by the corresponding sequence clipPath, as
		// well as another polygon (same points as clipPath polygon) with no fill to
		// draw the border.
		bar.selectAll(".sgroup")
			.data(data.order, s => s)
			.join(
				enter => {
					enter = enter.append("g")
						.attr("class", "sgroup")
						.attr("transform", translate)
					let inner = enter.append("g").attr("clip-path", s => `url(#${s}-clip)`)
					inner.append("rect")
						.attr("class", "seq-bg")
						.attr("fill", "white")
						.attr("width", s => x(data.synthases[s].sequence.length))
						.attr("height", s => y.bandwidth())
					inner.append("g")
						.selectAll("rect")
						.data(s => data.synthases[s].domains)
						.join(
							enter => enter.append("rect")
								.attr("x", d => x(d.start))
								.attr("class", d => `${d.type} seq-doms`)
								.attr("fill", d => colors(d.type))
								.attr("width", d => x(d.end - d.start))
								.attr("height", y.bandwidth())
								.on("mouseenter", cellEnter)
								.on("mouseleave", tooltipLeave),
						)
					enter.append("polygon")
						.attr("fill", "none")
						.attr("stroke", "black")
						.attr("stroke-width", "thin")
						.attr("points", s => points[s])
					return enter
				},
				update => update.call(
					update => {
						update.selectAll(".seq-bg")
							.transition(t)
							.attr("width", s => x(data.synthases[s].sequence.length))
						update.selectAll(".seq-doms")
							.transition(t)
							.attr("x", d => x(d.start))
							.attr("width", d => x(d.end - d.start))
						update.selectAll("polygon")
							.transition(t)
							.attr("points", s => points[s])
						return update.transition(t).attr("transform", translate)
					}
				),
			)

		// Adjust x axis on bottom
		const axisBottom = d3.axisBottom(x).ticks(6)
		xAxis.transition(t)
			.attr("transform", `translate(0, ${constants.barHeight * data.order.length})`)
			.call(axisBottom)

		// Adjust y axis labels on LHS of plot and add remove behaviour
		const axisLeft = d3.axisLeft(y)
			.tickSize(0)
			.tickPadding(10)
		yAxis.transition(t)
			.call(axisLeft)
		yAxis.selectAll("path").remove()

		function removeSequence(query) {
			// Create new data object.
			// 1. Remove sequence from order property
			// 2. Remove sequence from any classification groups (& remove empty)
			// 3. Remove sequence from synthases property
			let synth = data.synthases[query]
			let newData = {
				...data,
				"order": data.order.filter(q => q != query),
				"types": {
					...removeProps(data.types, synth.classification),
					...Object.fromEntries(
						synth.classification
							.map(c => [c, data.types[c].filter(q => q != query)])
							.filter(c => c[1].length > 0)
					)
				},
				"synthases": { ...removeProps(data.synthases, [query]) },
				"groups": data.groups.map(group => group.map(g => g))
			}
			// Prune classification groups:
			// 1. Delete specific classification if no sequences left
			// 2. Delete group in data.groups if now empty
			for (let [index, group] of newData.groups.entries()) {
				for (let [j, g] of group.entries()) {
					if (synth.classification.includes(g.classification)) {
						if (!newData.types.hasOwnProperty(g.classification))
							newData.groups[index].splice(j, 1)
						if (newData.groups[index].length == 0)
							newData.groups.splice(index, 1)
					}
				}
			}
			update(newData)
		}

		yAxis.selectAll("text")
			.style("font-size", "12px")
			.style("text-anchor", "end")
			.on("click", removeSequence)

		function annotate(groups) {
			/* Adds annotations for each classification group.
			*/
			let node = d3.select(this)
			let offsets = []
			previous = 0

			groups.forEach(group => {
				// Check for empty groups
				let synthases = data.types[group.classification]
				if (synthases.length == 0) return

				// If we've moved laterally, delete the last offset value.
				if (offsets && group.depth === previous)
					offsets.pop()

				// Calculate all the points of the bracket. Note that startX includes the
				// sum of all previous computed offsets, to allow for arbitrary level of
				// annotation nesting.
				let startX = (
					d3.max(synthases.map(s => x(data.synthases[s].sequence.length))) + 10
					+ d3.sum(offsets)
				)
				let endX = startX + 10
				let yPos = synthases.map(s => y(s))
				let topY = d3.min(yPos)
				let botY = d3.max(yPos) + y.bandwidth()
				let midY = botY + (topY - botY) / 2

				// Add a new group, containing the bracket path and the classification
				// text.
				let part = node.append("g")
				part.append("polyline")
					.attr("fill", "none")
					.attr("points", `${startX},${topY} ${endX},${topY} ${endX},${botY} ${startX},${botY}`)
					.attr("stroke", "black")
					.attr("stroke-width", "thin")
				part.append("text")
					.text(group.classification)
					.style("font-size", "12px")
					.attr("text-anchor", "start")
					.attr("x", endX + 10)
					.attr("y", midY)
					.attr("dy", "0.3em")

				// Add the total width of the added part to the array of previous offsets,
				// then set the depth.
				let offset = part.node().getBBox().width + 6
				offsets.push(offset)
				previous = group.depth
			})
		}

		// Add annotations
		annotations.selectAll("g")
			.data(data.groups, d => d.map(c => c.classification).join())
			.join("g")
			.each(annotate)

		// Draw the legend.
		// 1. Create scale based on original domain list for ordering
		//		This ensures ordering consistent with original legend (i.e. groups
		//		deleted and list shifted, rather than complete rearrangement)
		// 2. Sort current unique domains based on order scale
		// 3. Create height scale based on unique domains
		// 4. Draw legend groups
		let legDomain = getUniqueDomainProps(data).sort((a, b) => order(a) > order(b))
		let legScale = d3.scaleBand()
			.domain(legDomain)
			.range([0, constants.cellHeight * legDomain.length])
			.paddingInner(constants.cellPadding)
		let legTranslate = d => `translate(0, ${legScale(d)})`

		gLegend.selectAll("g")
			.data(legDomain, d => d)
			.join(
				enter => {
					enter = enter.append("g")
						.attr("transform", legTranslate)
					enter.append("rect")
						.attr("fill", d => colors(d))
						.attr("class", d => d)
						.attr("width", constants.cellWidth)
						.attr("height", legScale.bandwidth())
						.attr("cursor", "pointer")
						.on("click", e => {
							// Bind fill change behaviour then simulate click event
							picker.on("change", () => {
								d3.selectAll(`.${e}`).attr("fill", picker.node().value)
							})
							picker.node().click()
						})
					enter.append("text")
						.text(d => d)
						.attr("x", constants.cellWidth)
						.attr("y", legScale.bandwidth() / 2)
						.attr("dx", ".4em")
						.attr("dy", ".4em")
						.attr("text-anchor", "start")
						.style("font-size", "12px")
					return enter
				},
				update => update.call(
					update => update.transition(t).attr("transform", legTranslate)
				)
			)

		// Bind plot width input elements. Do inside update() so they always
		// reflect the current data (with removals)
		const widthTextInput = d3.select("#txt-plot-width")
		const widthRangeInput = d3.select("#rng-plot-width")
		widthTextInput.on("change", function() {
			widthRangeInput.property("value", this.value)
			constants.scaleFactor = +this.value
			update(data)
		})
		widthRangeInput.on("input", function() {
			widthTextInput.property("value", this.value)
			constants.scaleFactor = +this.value
			update(data)
		})

		// Adjust position of misc chart elements after transition has finished
		setTimeout(() => {
			const barBbox = bar.node().getBBox()
			title.transition(100)
				.text(`Domain architecture of ${data.order.length} sequences`)
				.attr("x", barBbox.width / 2)
			xLabel.transition(100)
				.attr("x", barBbox.width / 2)
				.attr("y", constants.barHeight * data.order.length + 40)
			let chartBbox = chart.node().getBBox()
			let legendBbox = gLegend.node().getBBox()
			let legendX = chartBbox.x + chartBbox.width + 20
			let legendY = chartBbox.y + chartBbox.height / 2 - legendBbox.height / 2
			gLegend.transition(100)
				.attr("transform", `translate(${legendX}, ${legendY})`)
		}, 300)
	}

	// Draw the initial plot
	update(data)
}
