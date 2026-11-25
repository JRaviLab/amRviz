#' Network UI Tab
#' @param id Module id; must match the id passed to `networkServer()`.
#' @return A tabPanel for drug-feature network visualization
#' @keywords internal
networkUI <- function(id) {
  ns <- shiny::NS(id)
  help <- function(label, text) {
    shiny::tags$span(class = "network-label", label, " ",
      shiny::tags$span(class = "network-help", tabindex = "0", role = "img", `aria-label` = text,
        "?", shiny::tags$span(class = "network-help-text", role = "tooltip", text)))
  }
  shiny::tabPanel(
    title = "Network",
    icon = shiny::icon("circle-nodes"),
    shiny::fluidPage(
      shiny::tags$style(shiny::HTML(.network_slider_css())),
      shiny::h3("Drug-Feature Network",
        style = paste0(
          "font-weight: bold; margin-top: 15px; margin-bottom: 15px;"
        )
      ),
      shiny::p(paste(
        "Interactive force-directed graph linking drugs or drug classes",
        "(dark diamonds) to their top predictive features of each type.",
        "Links are orange where a feature points toward resistance and blue",
        "toward susceptibility. Drag nodes, zoom, and hover for labels."
      )),
      shiny::radioButtons(
        ns("network_control_layout"), help("Controls", "Place controls above the network or in a left sidebar."),
        choices = c("Top" = "top", "Left" = "left"),
        selected = "top", inline = TRUE
      ),
      shiny::div(
        id = ns("network_workspace"), class = "network-workspace",
        shiny::div(class = "network-controls",
      shiny::fluidRow(
        class = "network-primary-controls",
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::selectInput(
            inputId = ns("network_bug_id"),
            label = shiny::tags$label(help("Bug", "Choose the species to display."), style = "font-size: 15px;"),
            choices = character(0),
            multiple = FALSE,
            selectize = TRUE,
            width = "100%"
          )
        ),
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::sliderInput(
            inputId = ns("network_top_n"),
            label = shiny::tags$label(
              help("Top features per type, per drug", "Keep this many features of each type per drug or drug class."),
              style = "font-size: 15px;"
            ),
            min = 1, max = 30, value = 2
          )
        ),
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          # Tell people why they can't use the dyad toggle when dyads are already forced on
          shiny::tags$div(
            id = ns("network_include_dyads_hint"),
            title = "Switch to 'Drugs \u2192 features' to use this option.",
            shiny::checkboxInput(
              inputId = ns("network_include_dyads"),
              label = help("Show dyads", "Dyads represent matched protein-gene pairs. Top dyads reflect pairs that occur repeatedly as important features across model seeds, helping to focus on potentially associated mechanisms."),
              value = FALSE
            )
          )
        ),
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::checkboxInput(ns("network_weighted_edges"), help("Weight edges by importance", "Vary edge thickness by feature importance or dyad contribution."), TRUE),
          shiny::tags$details(
            id = ns("network_feature_types_panel"),
            title = "In dyad view, enable 'Show connected features' to change feature types.",
            class = "network-feature-toggle",
            shiny::tags$summary(help("Toggle feature types", "Select feature nodes. In dyad view, enable Show connected features first. Choices do not filter the drug\u2013dyad links. Types without model support are grayed out.")),
            shiny::checkboxGroupInput(
              ns("network_feature_types"), label = NULL,
              choices = c("AMRFinder" = "args", "COGs" = "cogs", "Pfam" = "domains",
                          "Genes" = "genes", "Proteins" = "proteins", "Defense/Cas" = "DefenseCas",
                          "Structural" = "struct"),
              selected = names(SCALE_COLORS)
            )
          )
        )
        # ,
        # column(
        #   width = 3,
        #   style = "padding: 10px;",
        #   checkboxInput(
        #     inputId = "network_include_cogs",
        #     label = "Show COGs",
        #     value = FALSE
        #   ),
        #   tags$span(
        #     "Adds COG nodes from annotations.",
        #     style = "font-size: 11px; color: #666;"
        #   )
        # )
      ),
      shiny::fluidRow(
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::radioButtons(
            inputId = ns("network_rescored"),
            label = shiny::tags$label(
              help("Advanced: Rescored features and dyads", "Use saved scores aggregated across model seeds. Provides additional options for filtering and visualization."),
              style = "font-size: 15px;"
            ),
            choices = c("No" = "no", "Yes" = "yes"),
            selected = "no",
            inline = TRUE
          ),
          shiny::uiOutput(ns("network_rescoring_note")),
          shiny::conditionalPanel(
            condition = "input.network_rescored == 'yes'", ns = ns,
            shiny::radioButtons(
              ns("network_view"), help("Network view", "Connect drugs to features, or directly to dyads. Dyad edge color shows the resistance/susceptibility association; edge thickness follows contribution."),
              choices = c("Drugs \u2192 features" = "features", "Drugs \u2192 dyads" = "dyads"),
              selected = "dyads"
            ),
            shiny::conditionalPanel(
              condition = "input.network_view == 'dyads'", ns = ns,
              shiny::checkboxInput(
                ns("network_dyad_features"), help("Show connected features", "Add the selected feature types through their dyads. Hover over a dyad to see its supporting features."), FALSE
              )
            )
          )
        ),
        shiny::conditionalPanel(
          condition = "input.network_rescored == 'yes'",
          ns = ns,
          # One slider per category. Model performance is grey, features keep
          # Shiny's blue and dyads are green, the colour of dyad nodes. The
          # other rescoring thresholds keep everything (see rescoringAll()).
          shiny::column(
            width = 3,
            class = "network-sliders-grey",
            style = "padding: 10px;",
            shiny::tags$h5(help("Model performance", "Keep features from models that exceed this performance filter."), style = "font-weight: bold;"),
            shiny::sliderInput(
              inputId = ns("network_resc_delta"),
              label = shiny::tags$label(
                help("MCC above shuffled", "Models must beat their matched shuffled-label scores by at least this MCC. MCC ranges from -1 to +1; 0 is a random classifier."),
                style = "font-size: 15px;"
              ),
              min = 0, max = 1.2, value = rescoringAll()$delta_mcc, step = 0.05
            )
          ),
          shiny::column(
            width = 3,
            style = "padding: 10px;",
            shiny::tags$h5(help("Features", "Filter features by their rank scores across models and seeds."), style = "font-weight: bold;"),
            shiny::sliderInput(
              inputId = ns("network_resc_quantile"),
              label = shiny::tags$label(
                help("Rank-score", paste("Rank score is a feature's importance rank within a model, from 0 (least) to 1 (most), as the median across seeds.", "e.g., 0.75 keeps the top quarter by rank score.")),
                style = "font-size: 15px;"
              ),
              min = 0, max = 1, value = rescoringAll()$quantile, step = 0.01
            ),
            shiny::checkboxInput(
              inputId = ns("network_resc_both"),
              label = help("Found in binary and counts", "Keep a feature only if it is ranked in both binary and count models for the same feature type."),
              value = rescoringAll()$both_subtypes
            )
          ),
          shiny::column(
            width = 3,
            class = "network-sliders-green",
            style = "padding: 10px;",
            shiny::tags$h5(help("Dyads", "Filter rescored protein\u2013gene pairs."), style = "font-weight: bold;"),
            shiny::sliderInput(
              inputId = ns("network_resc_dyad"),
              label = shiny::tags$label(
                help("Dyad rank score", "The 0 to 1 range of rank scores for dyads available after Top dyads per drug filtering (below). The actual rank cutoff appears underneath this slider. Applies in dyad view or when Show dyads is on."),
                style = "font-size: 15px;"
              ),
              min = 0, max = 1, value = rescoringAll()$dyad_rank_score,
              step = 0.01
            ),
            shiny::tags$small(shiny::textOutput(ns("network_dyad_cutoff"), inline = TRUE)),
            shiny::conditionalPanel(
              condition = "input.network_view == 'dyads'", ns = ns,
              shiny::numericInput(
                ns("network_top_dyads"), help("Top dyads per drug", "Keep this many dyads per drug or class, ranked by score, then contribution."),
                value = 10, min = 1, step = 1
              )
            )
          )
        )
      )),
      shiny::fluidRow(
        class = "network-graph",
        shiny::column(
          width = 12,
          shiny::actionButton(
            ns("network_reset"), "Reset network", class = "network-reset",
            title = "Reset network", `aria-label` = "Reset network",
            icon = shiny::icon("rotate-right")
          ),
          shinycssloaders::withSpinner(
            networkD3::forceNetworkOutput(
              ns("drug_feature_network"),
              height = "600px"
            )
          ),
          style = "padding: 10px; border: 1px solid lightgray;"
        )
      ))
    )
  )
}


#' CSS that colours the grey and green groups of sliders
#'
#' Shiny draws sliders with ionRangeSlider; the skin's bar, handle and value
#' bubbles are recoloured inside each group.
#'
#' @return A CSS string.
#' @keywords internal
#' @noRd
.network_slider_css <- function() {
  one <- function(group, bar, dark) {
    sprintf(paste0(
      ".%1$s .irs--shiny .irs-bar { background: %2$s; border-color: %2$s; }",
      ".%1$s .irs--shiny .irs-handle { border-color: %3$s; }",
      ".%1$s .irs--shiny .irs-handle:hover, .%1$s .irs--shiny .irs-handle.state_hover",
      " { background-color: #f2f2f2; }",
      ".%1$s .irs--shiny .irs-from, .%1$s .irs--shiny .irs-to,",
      " .%1$s .irs--shiny .irs-single { background-color: %3$s; }",
      ".%1$s .irs--shiny .irs-from:before, .%1$s .irs--shiny .irs-to:before,",
      " .%1$s .irs--shiny .irs-single:before { border-top-color: %3$s; }"
    ), group, bar, dark)
  }
  paste(
    ".network-label { position: relative; display: inline-block; }",
    ".network-help { display: inline-flex; align-items: center; justify-content: center; width: 14px; height: 14px; border: 1px solid #aaa; border-radius: 50%; color: #777; font: 10px Arial, sans-serif; cursor: help; vertical-align: middle; }",
    ".network-help-text { display: none; position: absolute; top: 100%; left: 0; z-index: 1000; width: 260px; max-width: 80vw; padding: 8px; border: 1px solid #ccc; border-radius: 4px; background: #fff; color: #333; font: normal 12px/1.4 Arial, sans-serif; text-align: left; white-space: normal; box-shadow: 0 2px 6px #0002; }",
    ".network-help:hover .network-help-text, .network-help:focus .network-help-text { display: block; }",
    ".network-feature-toggle .shiny-options-group { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 0 12px; }",
    ".network-feature-toggle .checkbox { margin: 4px 0; }",
    ".network-primary-controls .shiny-input-container { margin-bottom: 0; }",
    ".network-primary-controls .shiny-input-container > .checkbox { margin: 0; }",
    ".network-graph { min-width: 0; }",
    ".network-reset { position: absolute; top: 12px; right: 16px; z-index: 5; border: 1px solid #ddd; border-radius: 3px; background: #fffd; color: #777; opacity: .65; padding: 4px 7px; font-size: 12px; }",
    ".network-reset:hover, .network-reset:focus { opacity: 1; }",
    "@media (min-width: 992px) {",
    " .network-workspace.network-controls-left { display: grid; grid-template-columns: 320px minmax(0, 1fr); gap: 20px; align-items: start; }",
    " .network-controls-left .network-controls > .row { margin: 0; }",
    " .network-controls-left .network-controls .col-sm-3 { width: 100%; float: none; }",
    " .network-controls-left .network-graph { margin: 0; position: sticky; top: 12px; }",
    " .network-controls-left .forceNetwork { height: var(--network-height, calc(100dvh - 260px)) !important; }",
    "}",
    ".network-feature-toggle { margin-top: 20px; }",
    ".network-feature-toggle > summary { list-style: none; cursor: pointer; font-weight: 600; }",
    ".network-feature-toggle[open] > summary { margin-bottom: 12px; }",
    ".network-feature-toggle > summary::-webkit-details-marker { display: none; }",
    ".network-feature-toggle > summary::before { content: '\\25B6'; display: inline-block; margin-right: 6px; }",
    ".network-feature-toggle[open] > summary::before { content: '\\25BC'; }",
    one("network-sliders-grey", "#a8a8a8", "#7a7a7a"),
    one("network-sliders-green", "#8ccf4a", .NETWORK_DYAD_COLOR)
  )
}

# Gently follow the canvas proportions in the desktop sidebar view. Reuse the
# widget's simulation, and disconnect its old observer on filter redraws.
.NETWORK_LAYOUT_JS <- "
function(el, x, state) {
  var force = this;
  if (!force || !force.force || !window.ResizeObserver) return;
  if (el.networkLayoutCleanup) el.networkLayoutCleanup();
  var svg = d3.select(el).select('svg'), frame, previous;
  function fit() {
    var nodes = force.nodes().filter(function(n) { return isFinite(n.x) && isFinite(n.y); });
    if (!nodes.length || !el.clientWidth || !el.clientHeight) return;
    var x0 = d3.min(nodes, function(n) { return n.x; });
    var x1 = d3.max(nodes, function(n) { return n.x; });
    var y0 = d3.min(nodes, function(n) { return n.y; });
    var y1 = d3.max(nodes, function(n) { return n.y; });
    var width = el.clientWidth, height = el.clientHeight, pad = 40;
    var k = Math.min(1, (width - 2 * pad) / Math.max(1, x1 - x0),
      (height - 2 * pad) / Math.max(1, y1 - y0));
    var transform = d3.zoomIdentity.translate(
      width / 2 - k * (x0 + x1) / 2, height / 2 - k * (y0 + y1) / 2).scale(k);
    svg.property('__zoom', transform);
    svg.select('.zoom-layer').attr('transform', transform);
  }
  var panel = el.closest('.network-graph');
  var button = panel && panel.querySelector('.network-reset');
  if (button) button.onclick = fit;
  function layout() {
    if (!el.clientWidth || !el.clientHeight) return;
    var workspace = el.closest('.network-workspace');
    var left = workspace && workspace.classList.contains('network-controls-left') &&
      window.matchMedia('(min-width: 992px)').matches;
    if (left) {
      var inset = el.getBoundingClientRect().top - panel.getBoundingClientRect().top;
      var top = Math.max(12, workspace.getBoundingClientRect().top) + inset;
      el.style.setProperty('--network-height', Math.max(160, Math.floor(window.innerHeight - top - 24)) + 'px');
    } else {
      el.style.removeProperty('--network-height');
    }
    var width = el.clientWidth, height = el.clientHeight;
    var key = [left, width, height].join(':');
    if (key === previous) return;
    previous = key;
    // Size the SVG immediately, including on the first sidebar render.
    svg.attr('width', width).attr('height', height);
    force.force('center', d3.forceCenter(width / 2, height / 2));
    var hadForces = !!force.force('layout-x');
    var aspect = width / height;
    var horizontal = 0.03 * Math.min(1, 1 / (aspect * aspect));
    var vertical = 0.02 * Math.min(1, aspect * aspect);
    force.force('layout-x', left ? d3.forceX(width / 2).strength(horizontal) : null);
    force.force('layout-y', left ? d3.forceY(height / 2).strength(vertical) : null);
    if (left || hadForces) force.alpha(0.4).restart();
  }
  // Wait for the default-filter graph to settle before fitting a full reset.
  if (state && state.reset && state.ready && el.networkResetApplied !== state.reset) {
    force.on('end.network-reset', function() {
      fit();
      el.networkResetApplied = state.reset;
      force.on('end.network-reset', null);
    });
  }
  function schedule() {
    cancelAnimationFrame(frame);
    frame = requestAnimationFrame(layout);
  }
  var observer = new ResizeObserver(schedule);
  observer.observe(el);
  window.addEventListener('resize', schedule);
  window.addEventListener('scroll', schedule, {passive: true});
  el.networkLayoutCleanup = function() {
    observer.disconnect();
    force.on('end.network-reset', null);
    cancelAnimationFrame(frame);
    window.removeEventListener('resize', schedule);
    window.removeEventListener('scroll', schedule);
  };
  layout();
}
"
