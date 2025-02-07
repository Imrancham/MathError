let parsedData = null; // Global variable to store parsed JSON data

// Handle JSON file upload
document.getElementById("jsonFileInput").addEventListener("change", function (event) {
  const file = event.target.files[0];
  if (!file) {
    alert("Please select a valid JSON file.");
    return;
  }

  const reader = new FileReader();
  reader.onload = function (e) {
    try {
      parsedData = JSON.parse(e.target.result);
      if (!Array.isArray(parsedData)) {
        throw new Error("Invalid JSON format. Expected an array of objects.");
      }
      alert("JSON file uploaded successfully! Now specify the range or IDs and click Search.");
    } catch (error) {
      console.error("Error parsing JSON:", error);
      parsedData = null;
      alert("Error processing the JSON file. Ensure it contains a valid array of objects.");
    }
  };
  reader.readAsText(file);
});

// Handle Search button click
document.getElementById("searchButton").addEventListener("click", function () {
  const start = parseInt(document.getElementById("startRow").value, 10);
  const end = parseInt(document.getElementById("endRow").value, 10);
  const idList = document.getElementById("idList").value.trim();

  if (!parsedData || !parsedData.length) {
    alert("Please upload a JSON file first.");
    return;
  }

  if (idList) {
    console.log("ID List:", idList);
    const ids = idList.trim().split(",").map((id) => id.trim());
    const filteredData = parsedData.filter((row) => ids.includes(String(row.id)));
    if (filteredData.length === 0) {
      alert("No data available for the specified IDs.");
      return; 
    }
    generateTable(filteredData);
    return;
  }

  if (isNaN(start) || isNaN(end) || start < 1 || end < start) {
    alert("Please enter a valid range (Start and End) or a list of IDs.");
    return;
  }

  const rangeData = parsedData.slice(start - 1, end);
  generateTable(rangeData);
});

function generateTable(data) {
  const container = document.getElementById("table-container");
  container.innerHTML = "";

  if (!data.length) {
    container.innerHTML = "<p>No data available for the specified criteria.</p>";
    return;
  }

  const table = document.createElement("table");
  const thead = document.createElement("thead");
  const tbody = document.createElement("tbody");

  const headers = Object.keys(data[0]);
  const headerRow = document.createElement("tr");
  headers.forEach((header) => {
    const th = document.createElement("th");
    th.textContent = header;
    headerRow.appendChild(th);
  });
  thead.appendChild(headerRow);

  const indexHeader = document.createElement("th");
  indexHeader.textContent = "Index";
  headerRow.insertBefore(indexHeader, headerRow.firstChild);

  data.forEach((row, index) => {
    const tr = document.createElement("tr");

    const indexCell = document.createElement("td");
    indexCell.textContent = index + 1;
    tr.appendChild(indexCell);

    headers.forEach((header) => {
      const td = document.createElement("td");
      if (header === "label") {
        const select = document.createElement("select");
        for (let i = -1; i <= 60; i++) {
          const option = document.createElement("option");
          option.value = i;
          option.textContent = i;
          if (parseInt(row[header], 10) === i) {
            option.selected = true;
          }
          select.appendChild(option);
        }
        select.addEventListener("change", function () {
          row.label = select.value;
          console.log("Updated row data:", row);
        });
        td.appendChild(select);
      } else if (["assumptions", "t0", "t1"].includes(header)) {
        const latexSpan = document.createElement("span");
        latexSpan.className = "latex-cell";
        katex.render(row[header] || "", latexSpan, { throwOnError: false });
        td.appendChild(latexSpan);
      } else {
        td.textContent = row[header];
      }
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });

  table.appendChild(thead);
  table.appendChild(tbody);
  container.appendChild(table);
}

// Handle Update JSON button click
document.getElementById("updateButton").addEventListener("click", function () {
  if (!parsedData || !parsedData.length) {
    alert("Please upload a JSON file first.");
    return;
  }

  const tableRows = document.querySelectorAll("#table-container table tbody tr");
  tableRows.forEach((row) => {
    const select = row.querySelector("select");
    const idCell = row.querySelector("td:nth-child(2)");
    if (select && idCell) {
      const id = idCell.textContent.trim();
      const parsedDataItem = parsedData.find(item => item.id === id);
      if (parsedDataItem) {
        parsedDataItem.label = select.value;
      }
    }
  });

  console.log("Updated JSON:", parsedData);

  const updatedJson = JSON.stringify(parsedData, null, 2);
  const blob = new Blob([updatedJson], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  const originalFileName = document.getElementById("jsonFileInput").files[0].name;
  const updatedFileName = originalFileName.replace(/\.json$/, '') + 'updated.json';
  link.download = updatedFileName;
  link.click();
});

// Add filter input and button
const filterContainer = document.createElement("div");
filterContainer.classList.add("input-container");

const filterInput = document.createElement("input");
filterInput.type = "number";
filterInput.id = "filterLabel";
filterInput.placeholder = "Filter by label (0-50)";
filterInput.min = 0;
filterInput.max = 50;

const filterButton = document.createElement("button");
filterButton.id = "filterButton";
filterButton.textContent = "Filter";

filterContainer.appendChild(filterInput);
filterContainer.appendChild(filterButton);
document.body.insertBefore(filterContainer, document.getElementById("table-container"));

// Handle Filter button click
document.getElementById("filterButton").addEventListener("click", function () {
  const filterValue = parseInt(document.getElementById("filterLabel").value, 10);
  if (isNaN(filterValue) || filterValue < 0 || filterValue > 50) {
    alert("Please enter a valid label value between 0 and 50.");
    return;
  }
  const start = parseInt(document.getElementById("startRow").value, 10);
  const end = parseInt(document.getElementById("endRow").value, 10);
  const idList = document.getElementById("idList").value.trim();

  if (!parsedData || !parsedData.length) {
    alert("Please upload a JSON file first.");
    return;
  }

  if (idList) {
    const ids = idList.split(",").map((id) => id.trim());
    const rangeData = parsedData.filter((row) => ids.includes(row.id));
    const filteredData = rangeData.filter((row) => parseInt(row.label, 10) === filterValue);

    generateTable(filteredData);
    return;
  }

  if (isNaN(start) || isNaN(end) || start < 1 || end < start) {
    alert("Please enter a valid range (Start and End) or a list of IDs.");
    return;
  }
  const rangeData = parsedData.slice(start - 1, end);

  const filteredData = rangeData.filter((row) => parseInt(row.label, 10) === filterValue);

  generateTable(filteredData);
});
