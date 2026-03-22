function sumValues(record) {
  return Object.values(record || {}).reduce((acc, value) => acc + (Number(value) || 0), 0);
}

function buildMetricItems(metrics) {
  const items = [];
  const queryCountByTaxid = metrics?.queryCountByTaxid || {};
  const fetchCountByTaxid = metrics?.fetchCountByTaxid || {};
  const boldSpecimenCountByTaxon = metrics?.boldSpecimenCountByTaxon || {};
  const boldDownloadedByTaxon = metrics?.boldDownloadedByTaxon || {};
  const boldMatchedByTaxon = metrics?.boldMatchedByTaxon || {};

  const allTaxids = new Set([...Object.keys(queryCountByTaxid), ...Object.keys(fetchCountByTaxid)]);
  for (const taxid of Array.from(allTaxids).sort()) {
    if (queryCountByTaxid[taxid] != null) {
      items.push(`query count taxid=${taxid}: ${queryCountByTaxid[taxid]}`);
    }
    if (fetchCountByTaxid[taxid] != null) {
      const fetched = fetchCountByTaxid[taxid];
      const total = queryCountByTaxid[taxid];
      items.push(
        total != null
          ? `fetch progress taxid=${taxid}: ${Math.min(fetched, total)}/${total}`
          : `fetch progress taxid=${taxid}: ${fetched}`
      );
    }
  }

  const allBoldTaxa = new Set([
    ...Object.keys(boldSpecimenCountByTaxon),
    ...Object.keys(boldDownloadedByTaxon),
    ...Object.keys(boldMatchedByTaxon)
  ]);
  for (const taxon of Array.from(allBoldTaxa).sort()) {
    const parts = [];
    if (boldSpecimenCountByTaxon[taxon] != null) {
      parts.push(`specimens=${boldSpecimenCountByTaxon[taxon]}`);
    }
    if (boldDownloadedByTaxon[taxon] != null) {
      parts.push(`downloaded=${boldDownloadedByTaxon[taxon]}`);
    }
    if (boldMatchedByTaxon[taxon] != null) {
      parts.push(`matched=${boldMatchedByTaxon[taxon]}`);
    }
    if (parts.length) {
      items.push(`bold taxon=${taxon}: ${parts.join(" ")}`);
    }
  }

  if (metrics?.matchedRecords != null) items.push(`matched records: ${metrics.matchedRecords}`);
  if (metrics?.keptRecordsBeforePostPrep != null) {
    items.push(`kept before post_prep: ${metrics.keptRecordsBeforePostPrep}`);
  }
  if (metrics?.primerTrimRemoved != null) items.push(`primer_trim removed: ${metrics.primerTrimRemoved}`);
  if (metrics?.lengthFilterRemoved != null) items.push(`length_filter removed: ${metrics.lengthFilterRemoved}`);
  if (metrics?.duplicateGroups != null) items.push(`duplicate groups: ${metrics.duplicateGroups}`);
  if (metrics?.crossOrganismGroups != null) {
    items.push(`cross_organism_groups: ${metrics.crossOrganismGroups}`);
  }

  return items;
}

function buildProgressDetail(phase, percent, metrics) {
  const safePercent = typeof percent === "number" ? Math.max(0, Math.min(percent, 100)) : null;
  const parts = [];

  const queryCountByTaxid = metrics?.queryCountByTaxid || {};
  const fetchCountByTaxid = metrics?.fetchCountByTaxid || {};
  const boldSpecimenCountByTaxon = metrics?.boldSpecimenCountByTaxon || {};
  const boldDownloadedByTaxon = metrics?.boldDownloadedByTaxon || {};
  const boldMatchedByTaxon = metrics?.boldMatchedByTaxon || {};

  const total = sumValues(queryCountByTaxid);
  if (total > 0) {
    const fetched = Object.entries(fetchCountByTaxid).reduce((acc, [taxid, raw]) => {
      const done = Number(raw) || 0;
      const max = Number(queryCountByTaxid[taxid]) || done;
      return acc + Math.min(done, max);
    }, 0);
    const safeFetched = Math.min(fetched, total);
    parts.push(`fetch ${safeFetched}/${total} (${((safeFetched / total) * 100).toFixed(1)}%)`);
  }

  const boldDownloadedTotal = sumValues(boldDownloadedByTaxon);
  const boldMatchedTotal = sumValues(boldMatchedByTaxon);
  const boldSpecimenTotal = sumValues(boldSpecimenCountByTaxon);
  if (boldDownloadedTotal > 0 || boldSpecimenTotal > 0) {
    const boldTaxonCount = new Set([
      ...Object.keys(boldSpecimenCountByTaxon),
      ...Object.keys(boldDownloadedByTaxon),
      ...Object.keys(boldMatchedByTaxon)
    ]).size;
    let detail = `bold matched ${boldMatchedTotal}/${boldDownloadedTotal || 0}`;
    if (boldSpecimenTotal > 0) {
      detail += ` specimens ${boldSpecimenTotal}`;
    }
    if (boldTaxonCount > 0) {
      detail += ` taxa ${boldTaxonCount}`;
    }
    parts.push(detail);
  }

  if (safePercent != null && phase !== "Fetch/Parse") {
    parts.push(`total ${safePercent.toFixed(1)}%`);
  }

  if (parts.length) return parts.join(" | ");
  if (safePercent != null) return `${safePercent.toFixed(1)}%`;
  return "-";
}

export function createMonitorView({
  statusEl,
  phaseEl,
  progressEl,
  progressDetailEl,
  logsEl,
  metricsListEl,
  resultFilesEl,
  openPath
}) {
  function reset() {
    statusEl.textContent = "Idle";
    phaseEl.textContent = "-";
    progressEl.value = 0;
    progressDetailEl.textContent = "0.0%";
    logsEl.textContent = "";
    metricsListEl.innerHTML = "";
  }

  function appendLog(line) {
    const current = logsEl.textContent;
    const next = current.length > 50000 ? current.slice(-30000) : current;
    logsEl.textContent = `${next}${line}\n`;
    logsEl.scrollTop = logsEl.scrollHeight;
  }

  function renderMetrics(metrics) {
    metricsListEl.innerHTML = "";
    for (const text of buildMetricItems(metrics)) {
      const li = document.createElement("li");
      li.textContent = text;
      metricsListEl.appendChild(li);
    }
  }

  function renderProgressDetail(phase, percent, metrics) {
    progressDetailEl.textContent = buildProgressDetail(phase, percent, metrics);
  }

  function renderResultFiles(files) {
    resultFilesEl.innerHTML = "";
    for (const path of files) {
      const li = document.createElement("li");
      const row = document.createElement("div");
      row.className = "file-row";

      const text = document.createElement("code");
      text.textContent = path;

      const button = document.createElement("button");
      button.textContent = "Open";
      button.addEventListener("click", async () => {
        await openPath(path);
      });

      row.append(text, button);
      li.appendChild(row);
      resultFilesEl.appendChild(li);
    }
  }

  return {
    reset,
    appendLog,
    renderMetrics,
    renderProgressDetail,
    renderResultFiles
  };
}
