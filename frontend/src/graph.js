import cytoscape from "cytoscape";

const STYLES = [
  {
    selector: "node",
    style: {
      label: "data(label)",
      "font-size": 12,
      "font-family": "Avenir Next, Segoe UI, sans-serif",
      color: "#1c1915",
      "text-valign": "bottom",
      "text-margin-y": 8,
      "text-background-color": "#fbf8f4",
      "text-background-opacity": 0.92,
      "text-background-padding": "3px",
      "text-max-width": 120,
      "text-wrap": "wrap",
      width: 28,
      height: 28,
      "border-width": 2,
      "border-color": "#fbf8f4",
    },
  },
  {
    selector: 'node[kind = "disease"]',
    style: {
      "background-color": "#0f6e62",
      width: 54,
      height: 54,
      "font-size": 13,
      "font-weight": 600,
      "text-max-width": 148,
    },
  },
  {
    selector: 'node[kind = "known_drug"]',
    style: { "background-color": "#1e4d7b", shape: "round-rectangle", width: 34, height: 24 },
  },
  {
    selector: 'node[kind = "candidate_drug"]',
    style: { "background-color": "#b86a1a", shape: "round-rectangle", width: 34, height: 24 },
  },
  {
    selector: 'node[kind = "contraindication"]',
    style: { "background-color": "#8d3a3a", shape: "diamond", width: 26, height: 26 },
  },
  {
    selector: 'node[kind = "protein"]',
    style: { "background-color": "#3e6b4a", shape: "ellipse" },
  },
  {
    selector: 'node[kind = "disease_protein"]',
    style: { "background-color": "#1f6b62", shape: "round-rectangle", width: 30, height: 30 },
  },
  {
    selector: 'node[kind = "pathway"]',
    style: { "background-color": "#5b4d78", shape: "hexagon", width: 30, height: 30 },
  },
  {
    selector: 'node[kind = "column"]',
    style: {
      "background-opacity": 0,
      "border-width": 0,
      width: 1,
      height: 1,
      "font-size": 13,
      "font-weight": 600,
      color: "#5c5348",
      "text-valign": "center",
      "text-halign": "center",
      "text-max-width": 150,
    },
  },
  {
    selector: "edge",
    style: {
      width: 1.6,
      "line-color": "#8d8478",
      "target-arrow-color": "#8d8478",
      "target-arrow-shape": "triangle",
      "curve-style": "straight",
      "text-margin-y": -10,
      label: "data(label)",
      "font-size": 10,
      color: "#4e473f",
      "text-background-color": "#f7f4ef",
      "text-background-opacity": 0.92,
      "text-background-padding": "2px",
      "text-rotation": "autorotate",
    },
  },
  {
    selector: 'edge[kind = "target"]',
    style: { width: 2, "line-color": "#3e6b4a", "target-arrow-color": "#3e6b4a" },
  },
  {
    selector: 'edge[kind = "membership"]',
    style: { width: 2, "line-color": "#5b4d78", "target-arrow-color": "#5b4d78" },
  },
  {
    selector: 'edge[kind = "association"]',
    style: {
      width: 1.6,
      "line-color": "#0f6e62",
      "target-arrow-color": "#0f6e62",
      "curve-style": "taxi",
      "taxi-direction": "downward",
      "taxi-turn": 160,
      label: "",
    },
  },
  {
    selector: 'edge[kind = "recorded"], edge[kind = "predicted"], edge[kind = "contraindication"]',
    style: {
      "curve-style": "taxi",
      "taxi-direction": "upward",
      "taxi-turn": 70,
      label: "",
    },
  },
  {
    selector: "edge.is-selected",
    style: { label: "data(label)", width: 3, "line-color": "#1c1915", "target-arrow-color": "#1c1915" },
  },
  {
    selector: 'edge[kind = "recorded"]',
    style: { width: 2.2, "line-color": "#1e4d7b", "target-arrow-color": "#1e4d7b" },
  },
  {
    selector: 'edge[kind = "predicted"]',
    style: {
      width: 2,
      "line-style": "dashed",
      "line-color": "#b86a1a",
      "target-arrow-color": "#b86a1a",
    },
  },
  {
    selector: 'edge[kind = "contraindication"]',
    style: { "line-color": "#8d3a3a", "target-arrow-color": "#8d3a3a", "line-style": "dotted" },
  },
  {
    selector: "node.is-selected",
    style: { "border-width": 4, "border-color": "#1c1915" },
  },
  {
    selector: "edge.is-selected",
    style: { width: 3, "line-color": "#1c1915", "target-arrow-color": "#1c1915" },
  },
  {
    selector: ".is-dim",
    style: { opacity: 0.12 },
  },
];

export function mountGraph(container, elements, onNode, onEdge) {
  const cy = cytoscape({
    container,
    elements,
    style: STYLES,
    layout: { name: "preset", fit: true, padding: 48 },
    wheelSensitivity: 0.18,
    minZoom: 0.25,
    maxZoom: 2.4,
  });
  cy.on("tap", "node", (event) => onNode(event.target));
  cy.on("tap", "edge", (event) => onEdge(event.target));
  cy.fit(undefined, 48);
  return cy;
}
