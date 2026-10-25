import { useEffect, useRef, useState } from "react";
import { mountGraph } from "./graph";

const KIND_LABEL = {
  disease: "Disease",
  known_drug: "Recorded treatment",
  candidate_drug: "Ranked candidate",
  protein: "Target protein",
  disease_protein: "Disease-associated protein",
  pathway: "Pathway",
  contraindication: "Contraindication",
  column: "Column",
};

function kindLabel(kind) {
  return KIND_LABEL[kind] || kind;
}

function percent(value) {
  if (value === undefined || value === null || Number.isNaN(value)) {
    return "—";
  }
  return `${(value * 100).toFixed(1)}%`;
}

async function getJson(url) {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`Request failed (${response.status})`);
  }
  return response.json();
}

export function App() {
  const stageRef = useRef(null);
  const cyRef = useRef(null);
  const openRef = useRef(() => {});
  const [metrics, setMetrics] = useState(null);
  const [diseases, setDiseases] = useState([]);
  const [query, setQuery] = useState("");
  const [view, setView] = useState(null);
  const [selected, setSelected] = useState(null);
  const [contraindications, setContraindications] = useState(false);
  const [showIndirect, setShowIndirect] = useState(false);
  const [needle, setNeedle] = useState("");
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");

  async function loadView(params, flags = {}) {
    setLoading(true);
    setError("");
    const contraFlag = flags.contraindications ?? contraindications;
    const indirectFlag = flags.indirect ?? showIndirect;
    const search = new URLSearchParams();
    if (params.q) search.set("q", params.q);
    if (params.disease_id) search.set("disease_id", params.disease_id);
    if (contraFlag) search.set("contraindications", "true");
    if (indirectFlag) search.set("indirect", "true");
    try {
      const payload = await getJson(`/api/view?${search.toString()}`);
      setView(payload);
      const diseaseNode = payload.graph?.nodes?.find((node) => node.data.kind === "disease");
      setSelected(diseaseNode ? { type: "node", id: diseaseNode.data.id } : null);
    } catch (cause) {
      setError(cause.message || "The graph could not be loaded.");
    } finally {
      setLoading(false);
    }
  }

  openRef.current = (diseaseId) => loadView({ disease_id: diseaseId });

  useEffect(() => {
    let cancelled = false;
    async function boot() {
      try {
        const [metricPayload, diseasePayload] = await Promise.all([
          getJson("/api/metrics"),
          getJson("/api/diseases"),
        ]);
        if (cancelled) return;
        setMetrics(metricPayload);
        setDiseases(diseasePayload.diseases);
        await loadView({});
      } catch (cause) {
        if (!cancelled) setError(cause.message || "The API is not responding.");
        setLoading(false);
      }
    }
    boot();
    return () => {
      cancelled = true;
    };
    // The first view is the disease overview.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  useEffect(() => {
    if (!stageRef.current || !view?.graph?.nodes?.length) {
      if (cyRef.current) {
        cyRef.current.destroy();
        cyRef.current = null;
      }
      return undefined;
    }
    const cy = mountGraph(
      stageRef.current,
      [...view.graph.nodes, ...view.graph.edges],
      (node) => {
        if (node.data("kind") === "column") return;
        if (view.mode === "overview" && node.data("kind") === "disease") {
          openRef.current(node.id());
          return;
        }
        setSelected({ type: "node", id: node.id() });
      },
      (edge) => setSelected({ type: "edge", id: edge.id() })
    );
    cyRef.current = cy;
    return () => {
      cy.destroy();
      cyRef.current = null;
    };
  }, [view]);

  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) return;
    cy.elements().removeClass("is-dim is-selected");
    if (selected) {
      const chosen = cy.getElementById(selected.id);
      if (chosen.nonempty()) chosen.addClass("is-selected");
    }
    const query = needle.trim().toLowerCase();
    if (query) {
      cy.nodes().forEach((node) => {
        if (node.data("kind") === "column") return;
        const label = `${node.data("full_label") || ""} ${node.id()}`.toLowerCase();
        if (!label.includes(query)) node.addClass("is-dim");
      });
      cy.edges().forEach((edge) => {
        if (edge.source().hasClass("is-dim") || edge.target().hasClass("is-dim")) {
          edge.addClass("is-dim");
        }
      });
      return;
    }
    const pathIds = selected?.type === "node"
      ? view?.graph?.nodes?.find((node) => node.data.id === selected.id)?.data.path_ids
      : null;
    if (!pathIds?.length) return;
    const keep = new Set(pathIds);
    cy.nodes().forEach((node) => {
      if (node.data("kind") === "column") return;
      if (!keep.has(node.id())) node.addClass("is-dim");
    });
    cy.edges().forEach((edge) => {
      if (!keep.has(edge.source().id()) || !keep.has(edge.target().id())) edge.addClass("is-dim");
    });
  }, [needle, selected, view]);

  function submitQuery(event) {
    event.preventDefault();
    loadView({ q: query });
  }

  function chooseNode(id) {
    setSelected({ type: "node", id });
    const cy = cyRef.current;
    const node = cy?.getElementById(id);
    if (!node || node.empty()) return;
    const pathIds = node.data("path_ids");
    if (pathIds?.length) {
      const focus = cy.nodes().filter((item) => pathIds.includes(item.id()));
      cy.animate({ fit: { eles: focus, padding: 64 } }, { duration: 240 });
      return;
    }
    cy.animate({ center: { eles: node } }, { duration: 240 });
  }

  const test = metrics?.test;
  const selectedNode =
    selected?.type === "node"
      ? view?.graph?.nodes?.find((node) => node.data.id === selected.id)?.data
      : null;
  const selectedEdge =
    selected?.type === "edge"
      ? view?.graph?.edges?.find((edge) => edge.data.id === selected.id)?.data
      : null;
  const activeDiseaseId = view?.disease_id;

  return (
    <div className="app">
      <header className="topbar">
        <div className="brand">
          <p className="eyebrow">Biomedical knowledge graph</p>
          <h1>ReMedy</h1>
        </div>
        <div className="metrics" aria-label="Held-out test metrics">
          <Metric
            label="Held-out accuracy"
            value={percent(test?.accuracy)}
            note={test ? `${test.correct} of ${test.total} pairs` : "Loading"}
          />
          <Metric label="AUPRC" value={test ? test.auprc.toFixed(3) : "—"} note="Primary rank metric" />
          <Metric label="AUROC" value={test ? test.auroc.toFixed(3) : "—"} note="Score separation" />
        </div>
        <p className="disclaimer">
          Research scores only. Not a treatment recommendation.
        </p>
      </header>

      <div className="workspace">
        <aside className="rail">
          <form onSubmit={submitQuery}>
            <label htmlFor="disease-query">Disease</label>
            <input
              id="disease-query"
              value={query}
              placeholder="Parkinson, epilepsy, ADHD…"
              onChange={(event) => setQuery(event.target.value)}
            />
            <button type="submit">Draw graph</button>
          </form>
          <ul className="disease-list">
            {diseases.map((disease) => (
              <li key={disease.id}>
                <button
                  type="button"
                  className={disease.id === activeDiseaseId ? "is-active" : ""}
                  onClick={() => {
                    setQuery(disease.label);
                    loadView({ disease_id: disease.id });
                  }}
                >
                  {disease.label}
                </button>
              </li>
            ))}
          </ul>
        </aside>

        <main className="stage-wrap">
          <div className="stage-head">
            <div>
              <h2>{view?.title || "Disease graph"}</h2>
              <p>{view?.summary || "Type a disease to draw its graph."}</p>
              {view?.mode === "focus" && view.omissions?.length > 0 && (
                <details className="omissions">
                  <summary>What this picture does not claim</summary>
                  <ul>
                    {view.omissions.map((item) => (
                      <li key={item}>{item}</li>
                    ))}
                  </ul>
                </details>
              )}
            </div>
            <div className="stage-actions">
              {view?.mode === "focus" && (
                <>
                  <input
                    className="node-filter"
                    value={needle}
                    placeholder="Find a drug or protein"
                    aria-label="Find a drug or protein"
                    onChange={(event) => setNeedle(event.target.value)}
                  />
                  <label className="check">
                    <input
                      type="checkbox"
                      checked={showIndirect}
                      onChange={(event) => {
                        const next = event.target.checked;
                        setShowIndirect(next);
                        if (activeDiseaseId) {
                          loadView(
                            { disease_id: activeDiseaseId },
                            { contraindications, indirect: next }
                          );
                        }
                      }}
                    />
                    Show indirect pathway bridges
                  </label>
                  <label className="check">
                    <input
                      type="checkbox"
                      checked={contraindications}
                      onChange={(event) => {
                        const next = event.target.checked;
                        setContraindications(next);
                        if (activeDiseaseId) {
                          loadView(
                            { disease_id: activeDiseaseId },
                            { contraindications: next, indirect: showIndirect }
                          );
                        }
                      }}
                    />
                    Show contraindications
                    {view.contraindication_total ? ` (${view.contraindication_total})` : ""}
                  </label>
                </>
              )}
              <button
                type="button"
                className="ghost"
                onClick={() => cyRef.current?.fit(undefined, 48)}
              >
                Fit graph
              </button>
            </div>
          </div>
          <div className="stage" ref={stageRef} role="img" aria-label="Disease knowledge graph" />
          {loading && <p className="overlay">Drawing the graph…</p>}
          {error && <p className="overlay error">{error}</p>}
          {view?.mode === "ambiguous" && (
            <div className="overlay card">
              <p>{view.summary}</p>
              {view.matches?.map((disease) => (
                <button
                  key={disease.id}
                  type="button"
                  onClick={() => loadView({ disease_id: disease.id })}
                >
                  {disease.label}
                </button>
              ))}
            </div>
          )}
          <ul className="legend" aria-label="Graph legend">
            <li><i className="swatch disease" /> Disease</li>
            <li><i className="swatch known" /> Recorded treatment</li>
            <li><i className="swatch candidate" /> Ranked candidate</li>
            <li><i className="swatch protein" /> Target protein</li>
            <li><i className="swatch disease-protein" /> Disease protein</li>
            <li><i className="swatch pathway" /> Pathway</li>
            <li><i className="swatch line target" /> Drug targets protein</li>
            <li><i className="swatch line membership" /> Protein in pathway</li>
            <li><i className="swatch line association" /> Disease to protein</li>
            <li><i className="swatch line recorded" /> Recorded use</li>
            <li><i className="swatch dashed" /> Rank score, not a treatment</li>
            <li><i className="swatch line safety" /> Contraindication</li>
          </ul>
        </main>

        <aside className="inspector">
          <Inspector node={selectedNode} edge={selectedEdge} view={view} onChoose={chooseNode} />
          {metrics && (
            <details className="ceiling">
              <summary>How the accuracy was measured</summary>
              <p>
                The cutoff is {metrics.threshold.toFixed(2)}, fixed before the test pairs were scored.
                Validation AUPRC chose the model. On these test scores, the best possible cutoff
                reaches {percent(metrics.test.oracle_accuracy)}. That ceiling is below 70%, so the
                screen reports the held-out result instead of a tuned figure.
              </p>
              <p>
                Biological similarity on the same test pairs: AUPRC{" "}
                {metrics.similarity_baseline_test.auprc.toFixed(3)}, accuracy{" "}
                {percent(metrics.similarity_baseline_test.accuracy)}.
              </p>
            </details>
          )}
        </aside>
      </div>
    </div>
  );
}

function Metric({ label, value, note }) {
  return (
    <div>
      <p className="metric-label">{label}</p>
      <p className="metric-value">{value}</p>
      <p className="metric-note">{note}</p>
    </div>
  );
}

function formatEvidence(item) {
  if (typeof item === "string") return item;
  return [item.source, item.activity, item.assay, item.publication, item.organism]
    .filter((value) => value && String(value).toLowerCase() !== "unknown")
    .join(" · ");
}

function EvidenceList({ items }) {
  if (!items?.length) return null;
  return (
    <ul className="evidence">
      {items.map((item, index) => (
        <li key={`${formatEvidence(item)}-${index}`}>{formatEvidence(item)}</li>
      ))}
    </ul>
  );
}

function Inspector({ node, edge, view, onChoose }) {
  if (edge) {
    return (
      <section>
        <p className="eyebrow">Relation</p>
        <h3>{edge.label}</h3>
        <p>
          {edge.source_type || "source"} {edge.source_id} → {edge.target_type || "target"} {edge.target_id}
        </p>
        {edge.note && <p>{edge.note}</p>}
        {edge.evidence_count > 0 && <p>{edge.evidence_count} evidence records on this edge</p>}
        {edge.evidence_count > 0 && !edge.evidence?.length && (
          <p>No evidence row was attached for this edge.</p>
        )}
        <EvidenceList items={edge.evidence} />
        {view?.mode === "focus" && <FocusLists view={view} onChoose={onChoose} />}
      </section>
    );
  }

  if (!node) {
    return (
      <section>
        <p className="eyebrow">Inspector</p>
        <h3>Choose a disease</h3>
        <p>The graph opens on the twelve diseases. Select one, or type a name.</p>
      </section>
    );
  }

  return (
    <section>
      <p className="eyebrow">{kindLabel(node.kind)}</p>
      <h3>{node.full_label}</h3>
      {node.score !== undefined && (
        <p className={node.predicted ? "score high" : "score"}>
          Score {Number(node.score).toFixed(3)}
          {node.predicted ? " · at or above 0.50" : " · below 0.50"}
        </p>
      )}
      {node.relation && <p>Recorded as {node.relation}</p>}
      {node.source_database && <p>Source database: {node.source_database}</p>}
      {node.detail && <p>{node.detail}</p>}
      {node.missing?.map((item) => (
        <p key={item}>{item}</p>
      ))}
      {node.paths?.map((path) => (
        <div key={path.title}>
          <h4>{path.title}</h4>
          {path.note && <p>{path.note}</p>}
          <ul className="evidence">
            {path.steps?.map((step) => (
              <li key={step}>{step}</li>
            ))}
          </ul>
          <EvidenceList items={path.evidence} />
        </div>
      ))}
      <NameList title="Shared proteins" values={node.shared_proteins} />
      <NameList title="Shared pathways" values={node.shared_pathways} />
      {view?.mode === "focus" && (
        <FocusLists view={view} onChoose={onChoose} />
      )}
    </section>
  );
}

function FocusLists({ view, onChoose }) {
  return (
    <>
      <h4>Ranked candidates</h4>
      <ol className="candidates">
        {view.candidates.map((candidate) => (
          <li key={candidate.drug_id}>
            <button type="button" onClick={() => onChoose(candidate.drug_id)}>
              <span>{candidate.drug_name}</span>
              <span className="score-figure">{candidate.score.toFixed(3)}</span>
            </button>
            <span className="bar" style={{ width: `${Math.round(candidate.score * 100)}%` }} />
          </li>
        ))}
      </ol>
      <h4>Recorded treatments on the graph</h4>
      <ul className="known-list">
        {view.known.map((drug) => (
          <li key={drug.drug_id}>
            <button type="button" onClick={() => onChoose(drug.drug_id)}>
              {drug.drug_name}
              <small>{drug.relation}</small>
            </button>
          </li>
        ))}
      </ul>
    </>
  );
}

function NameList({ title, values }) {
  if (!values?.length) return null;
  return (
    <>
      <h4>{title}</h4>
      <p>{values.join(" · ")}</p>
    </>
  );
}
