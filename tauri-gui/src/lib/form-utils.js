export const DEFAULTS = {
  outputPrefix: "MiFish",
  ncbiDb: "nucleotide",
  ncbiRettype: "gb",
  ncbiRetmode: "text",
  ncbiPerQuery: 100,
  defaultHeaderFormat: "{acc_id}|{organism}|{marker}|{label}|{type}|{loc}|{strand}",
  mifishHeaderFormat: "{db}|{acc_id}|{organism}"
};

export function parseIntOrNull(value) {
  const v = `${value}`.trim();
  if (!v) return null;
  const parsed = Number.parseInt(v, 10);
  return Number.isNaN(parsed) ? null : parsed;
}

export function parseFloatOrNull(value) {
  const v = `${value}`.trim();
  if (!v) return null;
  const parsed = Number.parseFloat(v);
  return Number.isNaN(parsed) ? null : parsed;
}

export function parseCommaSeparatedList(value) {
  return `${value}`
    .split(",")
    .map((item) => item.trim())
    .filter(Boolean);
}

export function readCheckedValues(elements) {
  return elements.filter((el) => el.checked).map((el) => el.value);
}

export function setCheckedValues(elements, values) {
  const wanted = new Set(Array.isArray(values) ? values : []);
  elements.forEach((el) => {
    el.checked = wanted.has(el.value);
  });
}
