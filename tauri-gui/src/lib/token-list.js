export function parseDelimitedTokens(raw, pattern = /[\s,;]+/) {
  return `${raw}`
    .split(pattern)
    .map((value) => value.trim())
    .filter(Boolean);
}

export function mergeUniqueValues(existing, values) {
  const merged = new Set(existing);
  const items = Array.isArray(values) ? values : [values];
  for (const raw of items) {
    const value = `${raw}`.trim();
    if (value) merged.add(value);
  }
  return Array.from(merged);
}

export function renderTokenPills({
  container,
  countLabel,
  items,
  countSuffix,
  pillClass = "taxid-pill",
  onRemove
}) {
  container.innerHTML = "";
  container.classList.toggle("empty", items.length === 0);
  if (countLabel) {
    countLabel.textContent = `${items.length} ${countSuffix}`;
  }

  for (const item of items) {
    const li = document.createElement("li");
    li.className = pillClass;

    const text = document.createElement("span");
    text.textContent = item;

    const remove = document.createElement("button");
    remove.type = "button";
    remove.className = "taxid-remove";
    remove.textContent = "×";
    remove.title = `remove ${item}`;
    remove.addEventListener("click", () => onRemove(item));

    li.append(text, remove);
    container.appendChild(li);
  }
}
