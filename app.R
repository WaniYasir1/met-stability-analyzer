library(shiny)
library(dplyr)
library(ggplot2)
library(DT)

# ── Dummy MMET dataset ────────────────────────────────────────────────────────
make_dummy <- function() {
  set.seed(42)
  envs <- paste0("E", 1:5)
  gens <- paste0("G", 1:8)
  reps <- paste0("R", 1:3)
  expand.grid(ENV = envs, GEN = gens, REP = reps) %>%
    mutate(
      ENV = as.factor(ENV),
      GEN = as.factor(GEN),
      REP = as.factor(REP),
      YIELD = round(
        30 +
          as.numeric(ENV) * 3 +
          as.numeric(GEN) * 2 +
          rnorm(n(), 0, 3), 2
      ),
      BIOMASS = round(
        100 +
          as.numeric(ENV) * 8 +
          as.numeric(GEN) * 5 +
          rnorm(n(), 0, 8), 2
      ),
      PROTEIN = round(
        12 +
          as.numeric(GEN) * 0.5 -
          as.numeric(ENV) * 0.3 +
          rnorm(n(), 0, 1), 2
      )
    )
}

dummy_data <- make_dummy()

# ── Analysis helpers ───────────────────────────────────────────────────────────

get_traits <- function(df) {
  fixed <- c("ENV","GEN","REP")
  setdiff(names(df)[sapply(df, is.numeric)], fixed)
}

run_cd <- function(df, traits) {
  n_rep <- length(unique(df$REP))
  bind_rows(lapply(traits, function(tr) {
    form <- as.formula(paste(tr, "~ ENV + GEN + REP"))
    aov_m <- aov(form, data = df)
    at <- anova(aov_m)
    MSE <- at["Residuals", "Mean Sq"]
    dfe <- at["Residuals", "Df"]
    t_v <- qt(0.975, dfe)
    SEd <- sqrt(2 * MSE / n_rep)
    CD  <- SEd * t_v
    data.frame(Trait=tr, MSE=round(MSE,4), df_error=dfe,
               t_val=round(t_v,4), SEd=round(SEd,4),
               CD_5pct=round(CD,4))
  }))
}

add_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*",
                              ifelse(p < 0.1, ".", " ")))))
}

run_gxe_anova <- function(df, trait) {
  form  <- as.formula(paste(trait, "~ GEN * ENV + REP"))
  model <- aov(form, data = df)
  at    <- summary(model)[[1]]
  data.frame(
    Source  = rownames(at),
    Df      = at[,"Df"],
    SS      = round(at[,"Sum Sq"], 3),
    MS      = round(at[,"Mean Sq"], 3),
    F_value = round(at[,"F value"], 3),
    P_value = signif(at[,"Pr(>F)"], 3),
    Signif  = add_stars(at[,"Pr(>F)"])
  )
}

run_env_anova <- function(df, trait) {
  form  <- as.formula(paste(trait, "~ ENV"))
  model <- aov(form, data = df)
  at    <- summary(model)[[1]]
  data.frame(
    Source  = rownames(at),
    Df      = at[,"Df"],
    SS      = round(at[,"Sum Sq"], 3),
    MS      = round(at[,"Mean Sq"], 3),
    F_value = round(at[,"F value"], 3),
    P_value = signif(at[,"Pr(>F)"], 3),
    Signif  = add_stars(at[,"Pr(>F)"])
  )
}

eberhart_stability <- function(df, trait) {
  means_ge <- df %>%
    group_by(GEN, ENV) %>%
    summarise(mean_val = mean(.data[[trait]], na.rm=TRUE), .groups="drop")
  
  grand_mean <- mean(means_ge$mean_val)
  env_idx    <- means_ge %>%
    group_by(ENV) %>%
    summarise(env_mean = mean(mean_val)) %>%
    mutate(Ij = env_mean - grand_mean)
  
  data_reg <- means_ge %>%
    left_join(env_idx, by = "ENV")
  
  stability <- data_reg %>%
    group_by(GEN) %>%
    summarise(
      Mean    = round(mean(mean_val), 3),
      bi      = round(coef(lm(mean_val ~ Ij))[2], 4),
      s2d     = {
        m   <- lm(mean_val ~ Ij)
        MSe <- summary(m)$sigma^2
        round(MSe, 4)
      },
      .groups = "drop"
    ) %>%
    mutate(
      Stability_Class = case_when(
        abs(bi - 1) < 0.2 & s2d < 1 ~ "✅ Stable",
        bi > 1.2                     ~ "📈 Responsive",
        bi < 0.8                     ~ "📉 Below avg env",
        TRUE                          ~ "⚠️ Unstable"
      )
    )
  list(stability = stability, env_idx = env_idx, data_reg = data_reg)
}

# ── UI ─────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(
    tags$link(
      href = paste0(
        "https://fonts.googleapis.com/css2?",
        "family=Space+Grotesk:wght@300;400;500;600;700&",
        "family=JetBrains+Mono:wght@400;500&",
        "family=Playfair+Display:wght@700;800&display=swap"
      ),
      rel = "stylesheet"
    ),
    tags$style(HTML("
      :root {
        --bg:       #0d1117;
        --surface:  #161b22;
        --surface2: #21262d;
        --border:   #30363d;
        --gold:     #d4a017;
        --gold-dim: #8a6800;
        --green:    #3fb950;
        --red:      #f85149;
        --blue:     #58a6ff;
        --text:     #e6edf3;
        --muted:    #8b949e;
        --font:     'Space Grotesk', sans-serif;
        --mono:     'JetBrains Mono', monospace;
        --serif:    'Playfair Display', serif;
      }
      *, *::before, *::after { box-sizing: border-box; }
      body {
        background: var(--bg);
        color: var(--text);
        font-family: var(--font);
        margin: 0; padding: 0;
        min-height: 100vh;
      }
      /* scrollbar */
      ::-webkit-scrollbar { width: 6px; height: 6px; }
      ::-webkit-scrollbar-track { background: var(--bg); }
      ::-webkit-scrollbar-thumb { background: var(--border); border-radius: 3px; }

      /* ── Header ── */
      .app-header {
        background: linear-gradient(135deg, var(--surface) 0%, #1a1f2e 100%);
        border-bottom: 1px solid var(--border);
        padding: 0 40px;
        height: 64px;
        display: flex; align-items: center; gap: 20px;
        position: sticky; top: 0; z-index: 200;
        box-shadow: 0 4px 30px rgba(0,0,0,0.5);
      }
      .header-logo {
        width: 38px; height: 38px;
        background: linear-gradient(135deg, var(--gold), #f0c040);
        border-radius: 9px;
        display: flex; align-items: center; justify-content: center;
        font-size: 18px; flex-shrink: 0;
        box-shadow: 0 0 16px rgba(212,160,23,0.4);
      }
      .header-title {
        font-family: var(--serif);
        font-size: 20px; color: var(--text); font-weight: 700;
        letter-spacing: -0.3px;
      }
      .header-sub {
        font-size: 11.5px; color: var(--muted);
        font-family: var(--mono); margin-left: auto;
      }
      .header-badge {
        background: rgba(212,160,23,0.15);
        border: 1px solid var(--gold-dim);
        color: var(--gold);
        font-family: var(--mono);
        font-size: 10px; padding: 3px 10px;
        border-radius: 20px; letter-spacing: 0.5px;
      }

      /* ── Layout ── */
      .main-wrap { display: flex; min-height: calc(100vh - 64px); }

      /* ── Sidebar ── */
      .sidebar {
        width: 240px; min-width: 240px;
        background: var(--surface);
        border-right: 1px solid var(--border);
        padding: 20px 0;
        position: sticky; top: 64px;
        height: calc(100vh - 64px);
        overflow-y: auto;
      }
      .sb-section {
        padding: 14px 18px 6px;
        font-family: var(--mono);
        font-size: 9.5px; color: var(--gold-dim);
        letter-spacing: 2px; text-transform: uppercase;
        border-top: 1px solid var(--border);
        margin-top: 6px;
      }
      .sb-section:first-child { border-top: none; margin-top: 0; }
      .nav-item {
        display: flex; align-items: center; gap: 10px;
        width: 100%; padding: 9px 18px;
        background: none; border: none;
        border-left: 3px solid transparent;
        color: var(--muted);
        font-family: var(--font); font-size: 13px; font-weight: 400;
        text-align: left; cursor: pointer;
        transition: all 0.15s ease;
      }
      .nav-item:hover {
        color: var(--text); background: rgba(255,255,255,0.04);
        border-left-color: var(--border);
      }
      .nav-item.active {
        color: var(--gold); background: rgba(212,160,23,0.08);
        border-left-color: var(--gold); font-weight: 500;
      }
      .nav-item .ni-icon {
        width: 26px; height: 26px; border-radius: 7px;
        background: rgba(255,255,255,0.06);
        display: flex; align-items: center; justify-content: center;
        font-size: 13px; flex-shrink: 0;
      }

      /* ── Content ── */
      .content { flex: 1; padding: 32px 36px; max-width: 1080px; overflow-x: hidden; }

      /* ── Page header ── */
      .pg-title {
        font-family: var(--serif);
        font-size: 28px; font-weight: 800;
        color: var(--text); margin-bottom: 4px;
        letter-spacing: -0.5px;
      }
      .pg-sub { font-size: 13.5px; color: var(--muted); margin-bottom: 28px; }

      /* ── Cards ── */
      .card {
        background: var(--surface);
        border: 1px solid var(--border);
        border-radius: 12px;
        margin-bottom: 20px;
        overflow: hidden;
        transition: box-shadow 0.2s;
      }
      .card:hover { box-shadow: 0 0 0 1px var(--border), 0 8px 32px rgba(0,0,0,0.4); }
      .card-hd {
        padding: 14px 20px;
        border-bottom: 1px solid var(--border);
        display: flex; align-items: center; gap: 10px;
        background: rgba(255,255,255,0.02);
      }
      .card-hd-title {
        font-weight: 600; font-size: 14px; color: var(--text);
        flex: 1;
      }
      .chip {
        font-family: var(--mono); font-size: 10px;
        padding: 2px 8px; border-radius: 20px; font-weight: 500;
      }
      .chip-gold   { background: rgba(212,160,23,0.15); color: var(--gold); border: 1px solid var(--gold-dim); }
      .chip-green  { background: rgba(63,185,80,0.12); color: var(--green); border: 1px solid #2ea043; }
      .chip-blue   { background: rgba(88,166,255,0.12); color: var(--blue); border: 1px solid #1f6feb; }
      .chip-red    { background: rgba(248,81,73,0.12); color: var(--red); border: 1px solid #da3633; }
      .card-bd { padding: 20px; }

      /* ── Info box ── */
      .info-box {
        background: rgba(212,160,23,0.06);
        border: 1px solid rgba(212,160,23,0.2);
        border-left: 4px solid var(--gold);
        border-radius: 8px;
        padding: 14px 18px;
        margin-bottom: 18px;
        font-size: 13px;
      }
      .info-box h4 {
        font-family: var(--mono); font-size: 11px; font-weight: 500;
        color: var(--gold); letter-spacing: 1px;
        text-transform: uppercase; margin-bottom: 8px;
      }
      .info-box ul { padding-left: 18px; color: #b0bec5; line-height: 1.9; }
      .info-box code {
        font-family: var(--mono); background: rgba(212,160,23,0.12);
        color: var(--gold); padding: 1px 5px; border-radius: 4px; font-size: 11.5px;
      }

      /* ── Interpretation ── */
      .interp-box {
        background: rgba(63,185,80,0.06);
        border: 1px solid rgba(63,185,80,0.2);
        border-left: 4px solid var(--green);
        border-radius: 8px;
        padding: 14px 18px; margin-top: 16px;
        font-size: 13px;
      }
      .interp-box h4 {
        font-family: var(--mono); font-size: 11px; font-weight: 500;
        color: var(--green); letter-spacing: 1px;
        text-transform: uppercase; margin-bottom: 8px;
      }
      .interp-box ul { padding-left: 18px; color: #b0bec5; line-height: 1.9; }
      .interp-box p  { color: #b0bec5; line-height: 1.7; }

      /* ── Buttons ── */
      .btn-run {
        background: var(--gold) !important;
        color: #0d1117 !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 9px 24px !important;
        font-family: var(--font) !important;
        font-size: 13px !important; font-weight: 600 !important;
        cursor: pointer; transition: all 0.15s !important;
        box-shadow: 0 0 12px rgba(212,160,23,0.3) !important;
      }
      .btn-run:hover {
        background: #f0c040 !important;
        box-shadow: 0 0 20px rgba(212,160,23,0.5) !important;
        transform: translateY(-1px);
      }
      .btn-load {
        background: var(--surface2) !important;
        color: var(--text) !important;
        border: 1px solid var(--border) !important;
        border-radius: 8px !important;
        padding: 8px 18px !important;
        font-family: var(--font) !important;
        font-size: 13px !important;
        cursor: pointer; transition: all 0.15s !important;
      }
      .btn-load:hover {
        border-color: var(--gold) !important;
        color: var(--gold) !important;
      }
      .btn-dl {
        background: rgba(88,166,255,0.12) !important;
        color: var(--blue) !important;
        border: 1px solid #1f6feb !important;
        border-radius: 8px !important;
        padding: 7px 16px !important;
        font-family: var(--font) !important;
        font-size: 12px !important; font-weight: 500 !important;
        cursor: pointer; transition: all 0.15s !important;
        margin-right: 6px !important; margin-top: 4px !important;
      }
      .btn-dl:hover {
        background: rgba(88,166,255,0.22) !important;
        box-shadow: 0 0 10px rgba(88,166,255,0.2) !important;
      }

      /* ── Form ── */
      .form-control, .selectize-input {
        background: var(--surface2) !important;
        border: 1px solid var(--border) !important;
        color: var(--text) !important;
        border-radius: 7px !important;
        font-family: var(--font) !important;
        font-size: 13px !important;
      }
      .form-control:focus, .selectize-input.focus {
        border-color: var(--gold) !important;
        box-shadow: 0 0 0 3px rgba(212,160,23,0.15) !important;
      }
      label { color: var(--muted) !important; font-size: 12.5px !important; font-weight: 500 !important; }
      .selectize-dropdown {
        background: var(--surface2) !important;
        border: 1px solid var(--border) !important;
        border-radius: 7px !important;
      }
      .selectize-dropdown-content .option { color: var(--text) !important; font-size: 13px !important; }
      .selectize-dropdown-content .option:hover { background: rgba(212,160,23,0.1) !important; }

      /* ── Tabs ── */
      .nav-tabs {
        border-bottom: 1px solid var(--border) !important;
        margin-bottom: 0 !important;
      }
      .nav-tabs > li > a {
        background: none !important; border: none !important;
        border-bottom: 2px solid transparent !important;
        color: var(--muted) !important;
        font-family: var(--font) !important;
        font-size: 12.5px !important; font-weight: 500 !important;
        padding: 9px 16px !important; margin-bottom: -1px !important;
        border-radius: 0 !important; transition: all 0.15s !important;
      }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        background: none !important; border: none !important;
        border-bottom: 2px solid var(--gold) !important;
        color: var(--gold) !important;
      }
      .tab-content { padding-top: 18px; }

      /* ── Output ── */
      pre.shiny-text-output {
        background: var(--surface2) !important;
        color: #a8d5a2 !important;
        border: 1px solid var(--border) !important;
        border-radius: 8px !important;
        font-family: var(--mono) !important;
        font-size: 12px !important; line-height: 1.7 !important;
        padding: 16px !important;
        max-height: 420px; overflow-y: auto;
      }

      /* ── DT table ── */
      .dataTables_wrapper { color: var(--text) !important; }
      table.dataTable { color: var(--text) !important; background: transparent !important; }
      table.dataTable thead th {
        background: var(--surface2) !important;
        color: var(--muted) !important;
        border-bottom: 1px solid var(--border) !important;
        font-family: var(--mono) !important;
        font-size: 11px !important; font-weight: 500 !important;
        letter-spacing: 0.5px;
      }
      table.dataTable tbody tr { background: transparent !important; }
      table.dataTable tbody tr:hover td { background: rgba(212,160,23,0.06) !important; }
      table.dataTable tbody td {
        border-bottom: 1px solid rgba(48,54,61,0.5) !important;
        font-size: 12.5px !important;
      }
      .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate {
        color: var(--muted) !important; font-size: 12px !important;
      }
      .paginate_button { color: var(--muted) !important; }
      .paginate_button.current { background: rgba(212,160,23,0.15) !important; color: var(--gold) !important; border-radius: 5px !important; }

      /* ── File input ── */
      .btn-file { background: var(--surface2) !important; color: var(--text) !important; border-color: var(--border) !important; border-radius: 7px !important; }

      /* ── Stat cards ── */
      .stat-grid { display: grid; grid-template-columns: repeat(3,1fr); gap: 14px; margin-bottom: 20px; }
      .stat-card {
        background: var(--surface2);
        border: 1px solid var(--border);
        border-radius: 10px; padding: 16px;
      }
      .stat-card .sc-label { font-family: var(--mono); font-size: 10px; color: var(--muted); letter-spacing: 1px; text-transform: uppercase; margin-bottom: 6px; }
      .stat-card .sc-val  { font-size: 22px; font-weight: 700; color: var(--gold); font-family: var(--serif); }
      .stat-card .sc-sub  { font-size: 11px; color: var(--muted); margin-top: 2px; }

      /* ── Welcome grid ── */
      .welcome-grid { display: grid; grid-template-columns: repeat(auto-fill,minmax(190px,1fr)); gap: 14px; margin-top: 22px; }
      .wc-card {
        background: var(--surface);
        border: 1px solid var(--border);
        border-radius: 10px; padding: 18px;
        transition: all 0.18s; cursor: pointer;
      }
      .wc-card:hover { border-color: var(--gold); box-shadow: 0 0 20px rgba(212,160,23,0.12); transform: translateY(-2px); }
      .wc-icon { font-size: 26px; margin-bottom: 10px; }
      .wc-title { font-weight: 600; font-size: 13.5px; color: var(--text); margin-bottom: 4px; }
      .wc-desc  { font-size: 12px; color: var(--muted); line-height: 1.5; }

      hr.div { border: none; border-top: 1px solid var(--border); margin: 14px 0; }
      .row { margin: 0 !important; }
      .col-md-4, .col-md-8, .col-md-6, .col-md-12 { padding: 0 8px !important; }
    "))
  ),
  
  # Header
  div(class="app-header",
      div(class="header-logo","🌾"),
      div(class="header-title","MET Stability Analyzer"),
      span(class="header-badge","agricolae"),
      div(class="header-sub","Multi-Environment Trial • GxE Analysis Suite")
  ),
  
  div(class="main-wrap",
      
      # ── Sidebar ──────────────────────────────────────────────────────────────
      div(class="sidebar",
          div(class="sb-section","Navigation"),
          tags$button(class="nav-item active", id="nav_home",  onclick="goPage('home')",   div(class="ni-icon","🏠"), "Home"),
          tags$button(class="nav-item",        id="nav_data",  onclick="goPage('data')",   div(class="ni-icon","📂"), "Data Import"),
          div(class="sb-section","Analysis"),
          tags$button(class="nav-item", id="nav_cd",    onclick="goPage('cd')",    div(class="ni-icon","📐"), "Critical Difference"),
          tags$button(class="nav-item", id="nav_gxe",   onclick="goPage('gxe')",   div(class="ni-icon","⚗️"), "G×E ANOVA"),
          tags$button(class="nav-item", id="nav_env",   onclick="goPage('env')",   div(class="ni-icon","🌍"), "Environment ANOVA"),
          tags$button(class="nav-item", id="nav_stab",  onclick="goPage('stab')",  div(class="ni-icon","📈"), "Stability (Eberhart)"),
          div(class="sb-section","Output"),
          tags$button(class="nav-item", id="nav_export",onclick="goPage('export')",div(class="ni-icon","💾"), "Export All Results")
      ),
      
      div(class="content",
          
          # ── HOME ──────────────────────────────────────────────────────────────
          div(id="pg_home",
              div(class="pg-title","Multi-Environment Trial Analyzer"),
              div(class="pg-sub","Complete pipeline: Critical Difference → G×E ANOVA → Environment ANOVA → Eberhart-Russell Stability"),
              div(class="welcome-grid",
                  div(class="wc-card", div(class="wc-icon","📂"), div(class="wc-title","Data Import"), div(class="wc-desc","Upload CSV or load built-in dummy MET data")),
                  div(class="wc-card", div(class="wc-icon","📐"), div(class="wc-title","Critical Difference"), div(class="wc-desc","MSE-based CD values for all traits")),
                  div(class="wc-card", div(class="wc-icon","⚗️"), div(class="wc-title","G×E ANOVA"), div(class="wc-desc","Genotype × Environment interaction analysis")),
                  div(class="wc-card", div(class="wc-icon","🌍"), div(class="wc-title","Env ANOVA"), div(class="wc-desc","Environment-only ANOVA for each trait")),
                  div(class="wc-card", div(class="wc-icon","📈"), div(class="wc-title","Stability"), div(class="wc-desc","Eberhart-Russell regression stability statistics")),
                  div(class="wc-card", div(class="wc-icon","💾"), div(class="wc-title","Export"), div(class="wc-desc","Download all results as TXT, CSV, and plots as PNG"))
              )
          ),
          
          # ── DATA IMPORT ───────────────────────────────────────────────────────
          div(id="pg_data", style="display:none",
              div(class="pg-title","Data Import"),
              div(class="pg-sub","Upload your MET dataset or use the built-in dummy data to explore"),
              div(class="info-box",
                  tags$h4("📋 Required Data Format"),
                  tags$ul(
                    tags$li("CSV file with columns: ", tags$code("ENV"), ", ", tags$code("GEN"), ", ", tags$code("REP"), ", and one or more numeric trait columns"),
                    tags$li(tags$code("ENV"), ": Environment names (e.g., E1, E2, E3...)"),
                    tags$li(tags$code("GEN"), ": Genotype names (e.g., G1, G2, G3...)"),
                    tags$li(tags$code("REP"), ": Replicate numbers (e.g., R1, R2, R3)"),
                    tags$li("Trait columns (e.g., ", tags$code("YIELD"), ", ", tags$code("BIOMASS"), ", ", tags$code("PROTEIN"), "): numeric, one column per trait"),
                    tags$li("All factor combinations (ENV × GEN × REP) must be balanced for best results")
                  )
              ),
              div(class="card",
                  div(class="card-hd", span(class="card-hd-title","Load Data"), span(class="chip chip-gold","Import")),
                  div(class="card-bd",
                      fluidRow(
                        column(5,
                               fileInput("upload","Upload CSV File", accept=".csv"),
                               div(style="margin-top:-10px; margin-bottom:14px;"),
                               actionButton("load_dummy","📥 Load Dummy Data", class="btn-load"),
                               p(style="font-size:11.5px; color:var(--muted); margin-top:8px;",
                                 "Dummy data: 8 genotypes × 5 environments × 3 replicates, traits: YIELD, BIOMASS, PROTEIN")
                        ),
                        column(7,
                               uiOutput("data_summary_ui"),
                               br(),
                               DTOutput("data_preview")
                        )
                      )
                  )
              )
          ),
          
          # ── CRITICAL DIFFERENCE ───────────────────────────────────────────────
          div(id="pg_cd", style="display:none",
              div(class="pg-title","Critical Difference (CD)"),
              div(class="pg-sub","Minimum difference between two treatment means to be considered significant at 5%"),
              div(class="info-box",
                  tags$h4("📋 How It Works"),
                  tags$ul(
                    tags$li("CD = t(α, df_error) × SEd   where   SEd = √(2 × MSE / r)"),
                    tags$li(tags$code("MSE"), ": Mean Square Error from ANOVA | ", tags$code("r"), ": Number of replicates"),
                    tags$li("If the difference between two genotype means > CD, they are significantly different"),
                    tags$li("Computed separately for each numeric trait in your dataset")
                  )
              ),
              div(class="card",
                  div(class="card-hd", span(class="card-hd-title","Run Critical Difference"), span(class="chip chip-gold","CD @ 5%")),
                  div(class="card-bd",
                      actionButton("run_cd","▶ Compute CD", class="btn-run"),
                      br(), br(),
                      tabsetPanel(
                        tabPanel("📊 CD Table", DTOutput("cd_table")),
                        tabPanel("📈 CD Plot",  plotOutput("cd_plot", height="360px"))
                      ),
                      uiOutput("cd_interp")
                  )
              )
          ),
          
          # ── G×E ANOVA ─────────────────────────────────────────────────────────
          div(id="pg_gxe", style="display:none",
              div(class="pg-title","G×E Interaction ANOVA"),
              div(class="pg-sub","Tests significance of Genotype, Environment, and their interaction effects"),
              div(class="info-box",
                  tags$h4("📋 Model"),
                  tags$ul(
                    tags$li("Model: ", tags$code("Trait ~ GEN * ENV + REP")),
                    tags$li("Significant GEN × ENV interaction means genotype ranking changes across environments"),
                    tags$li("Significance codes: *** p<0.001  ** p<0.01  * p<0.05  . p<0.1"),
                    tags$li("Run this for each trait; select trait from dropdown below")
                  )
              ),
              div(class="card",
                  div(class="card-hd", span(class="card-hd-title","G×E ANOVA"), span(class="chip chip-blue","Interaction")),
                  div(class="card-bd",
                      fluidRow(
                        column(4, uiOutput("gxe_trait_ui")),
                        column(8, br(), actionButton("run_gxe","▶ Run G×E ANOVA", class="btn-run"))
                      ),
                      br(),
                      tabsetPanel(
                        tabPanel("📊 ANOVA Table", DTOutput("gxe_table")),
                        tabPanel("📈 Interaction Plot", plotOutput("gxe_plot", height="380px"))
                      ),
                      uiOutput("gxe_interp")
                  )
              )
          ),
          
          # ── ENV ANOVA ─────────────────────────────────────────────────────────
          div(id="pg_env", style="display:none",
              div(class="pg-title","Environment-Only ANOVA"),
              div(class="pg-sub","Tests whether environments differ significantly for each trait"),
              div(class="info-box",
                  tags$h4("📋 Model"),
                  tags$ul(
                    tags$li("Simple one-way model: ", tags$code("Trait ~ ENV")),
                    tags$li("Significant result → environments create different yield potential"),
                    tags$li("Used to justify multi-environment testing and stability analysis"),
                    tags$li("All traits are analyzed automatically")
                  )
              ),
              div(class="card",
                  div(class="card-hd", span(class="card-hd-title","Environment ANOVA"), span(class="chip chip-green","Env Effect")),
                  div(class="card-bd",
                      actionButton("run_env","▶ Run Env ANOVA", class="btn-run"),
                      br(), br(),
                      uiOutput("env_trait_tabs"),
                      uiOutput("env_interp")
                  )
              )
          ),
          
          # ── STABILITY ─────────────────────────────────────────────────────────
          div(id="pg_stab", style="display:none",
              div(class="pg-title","Eberhart-Russell Stability Analysis"),
              div(class="pg-sub","Regression-based stability: bi (regression coefficient) and s²d (deviation from regression)"),
              div(class="info-box",
                  tags$h4("📋 Interpretation Guide"),
                  tags$ul(
                    tags$li(tags$code("Mean"), ": Average performance of genotype across all environments"),
                    tags$li(tags$code("bi ≈ 1"), ": Average stability — responds proportionally to environment quality"),
                    tags$li(tags$code("bi > 1"), ": Responsive — thrives in favorable environments; poor in unfavorable"),
                    tags$li(tags$code("bi < 1"), ": Stable — consistent in poor environments; unresponsive to good conditions"),
                    tags$li(tags$code("s²d ≈ 0"), ": High predictability — deviations from regression line are minimal"),
                    tags$li("✅ Ideal stable cultivar: High mean + bi ≈ 1 + s²d ≈ 0")
                  )
              ),
              div(class="card",
                  div(class="card-hd", span(class="card-hd-title","Stability Statistics"), span(class="chip chip-gold","Eberhart-Russell")),
                  div(class="card-bd",
                      fluidRow(
                        column(4, uiOutput("stab_trait_ui")),
                        column(8, br(), actionButton("run_stab","▶ Run Stability", class="btn-run"))
                      ),
                      br(),
                      tabsetPanel(
                        tabPanel("📊 Stability Table",   DTOutput("stab_table")),
                        tabPanel("📈 Stability Plot",     plotOutput("stab_plot",  height="400px")),
                        tabPanel("📈 Regression Plot",    plotOutput("stab_regplot",height="400px"))
                      ),
                      uiOutput("stab_interp")
                  )
              )
          ),
          
          # ── EXPORT ────────────────────────────────────────────────────────────
          div(id="pg_export", style="display:none",
              div(class="pg-title","Export All Results"),
              div(class="pg-sub","Download analysis results and plots in your preferred format"),
              div(class="card",
                  div(class="card-hd", span(class="card-hd-title","Download Options"), span(class="chip chip-blue","Export")),
                  div(class="card-bd",
                      p(style="color:var(--muted); font-size:13px; margin-bottom:18px;",
                        "Run the analyses first (CD, G×E ANOVA, Env ANOVA, Stability) before downloading."),
                      div(class="card", style="background:var(--surface2); border-color:var(--border);",
                          div(class="card-hd", span(class="card-hd-title","📄 Text Reports")),
                          div(class="card-bd",
                              downloadButton("dl_cd_csv",   "⬇ CD Results (.csv)",         class="btn-dl"),
                              downloadButton("dl_gxe_csv",  "⬇ G×E ANOVA (.csv)",          class="btn-dl"),
                              downloadButton("dl_env_csv",  "⬇ Env ANOVA (.csv)",           class="btn-dl"),
                              downloadButton("dl_stab_csv", "⬇ Stability Results (.csv)",   class="btn-dl"),
                              downloadButton("dl_all_txt",  "⬇ All Results (.txt)",          class="btn-dl")
                          )
                      ),
                      div(class="card", style="background:var(--surface2); border-color:var(--border);",
                          div(class="card-hd", span(class="card-hd-title","🖼 Plots")),
                          div(class="card-bd",
                              downloadButton("dl_cd_plot",   "⬇ CD Plot (.png)",              class="btn-dl"),
                              downloadButton("dl_stab_plot", "⬇ Stability Plot (.png)",       class="btn-dl"),
                              downloadButton("dl_reg_plot",  "⬇ Regression Plot (.png)",      class="btn-dl"),
                              downloadButton("dl_gxe_plot",  "⬇ G×E Interaction Plot (.png)", class="btn-dl"),
                              downloadButton("dl_data_csv",  "⬇ Download Data (.csv)",        class="btn-dl")
                          )
                      )
                  )
              )
          )
          
      ) # content
  ), # main-wrap
  
  tags$script(HTML("
    var allPages = ['home','data','cd','gxe','env','stab','export'];
    function goPage(p) {
      allPages.forEach(function(pg) {
        var el = document.getElementById('pg_'+pg);
        if(el) el.style.display = (pg===p) ? 'block' : 'none';
        var btn = document.getElementById('nav_'+pg);
        if(btn) btn.className = 'nav-item' + (pg===p ? ' active' : '');
      });
    }
  "))
)

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── Reactive data ──────────────────────────────────────────────────────────
  rv <- reactiveValues(df = NULL, traits = NULL,
                       cd = NULL, gxe = NULL, env_res = NULL,
                       stab = NULL, stab_env = NULL, stab_reg = NULL,
                       selected_gxe_trait = NULL,
                       selected_stab_trait = NULL)
  
  observeEvent(input$load_dummy, {
    rv$df     <- dummy_data
    rv$traits <- get_traits(dummy_data)
  })
  
  observeEvent(input$upload, {
    req(input$upload)
    df <- read.csv(input$upload$datapath, stringsAsFactors = FALSE)
    df <- df %>% mutate(across(any_of(c("ENV","GEN","REP")), as.factor))
    rv$df     <- df
    rv$traits <- get_traits(df)
  })
  
  # ── Data summary ────────────────────────────────────────────────────────────
  output$data_summary_ui <- renderUI({
    req(rv$df)
    df <- rv$df
    div(class="stat-grid",
        div(class="stat-card", div(class="sc-label","Environments"), div(class="sc-val", length(unique(df$ENV))), div(class="sc-sub","unique ENVs")),
        div(class="stat-card", div(class="sc-label","Genotypes"),    div(class="sc-val", length(unique(df$GEN))), div(class="sc-sub","unique GENs")),
        div(class="stat-card", div(class="sc-label","Traits"),       div(class="sc-val", length(rv$traits)),      div(class="sc-sub",paste(rv$traits, collapse=", ")))
    )
  })
  
  output$data_preview <- renderDT({
    req(rv$df)
    datatable(head(rv$df, 30), options=list(pageLength=10, scrollX=TRUE), rownames=FALSE)
  })
  
  # ── Trait selectors ─────────────────────────────────────────────────────────
  output$gxe_trait_ui <- renderUI({
    req(rv$traits)
    selectInput("gxe_trait","Select Trait:", choices=rv$traits, selected=rv$traits[1])
  })
  output$stab_trait_ui <- renderUI({
    req(rv$traits)
    selectInput("stab_trait","Select Trait:", choices=rv$traits, selected=rv$traits[1])
  })
  
  # ── CD ───────────────────────────────────────────────────────────────────────
  observeEvent(input$run_cd, {
    req(rv$df, rv$traits)
    rv$cd <- run_cd(rv$df, rv$traits)
  })
  
  output$cd_table <- renderDT({
    req(rv$cd)
    datatable(rv$cd, options=list(pageLength=15), rownames=FALSE) %>%
      formatRound(c("MSE","t_val","SEd","CD_5pct"), 4)
  })
  
  cd_plot_obj <- reactive({
    req(rv$cd)
    ggplot(rv$cd, aes(x=Trait, y=CD_5pct, fill=Trait)) +
      geom_bar(stat="identity", width=0.55, alpha=0.9) +
      geom_text(aes(label=round(CD_5pct,3)), vjust=-0.4, size=4, color="#e6edf3", fontface="bold") +
      scale_fill_manual(values=c("#d4a017","#3fb950","#58a6ff","#f85149","#bc8cff","#ff9f1c")) +
      labs(title="Critical Difference @ 5% by Trait", x="Trait", y="CD Value") +
      theme_minimal(base_size=13) +
      theme(
        plot.background=element_rect(fill="#161b22",color=NA),
        panel.background=element_rect(fill="#161b22",color=NA),
        panel.grid.major=element_line(color="#21262d"),
        panel.grid.minor=element_blank(),
        axis.text=element_text(color="#8b949e"),
        axis.title=element_text(color="#8b949e"),
        plot.title=element_text(color="#e6edf3",face="bold",size=14),
        legend.position="none"
      )
  })
  
  output$cd_plot <- renderPlot({ cd_plot_obj() })
  
  output$cd_interp <- renderUI({
    req(rv$cd)
    best <- rv$cd$Trait[which.min(rv$cd$CD_5pct)]
    div(class="interp-box",
        tags$h4("💡 Interpretation"),
        tags$ul(
          tags$li("CD (Critical Difference) is the minimum difference between two genotype means needed to declare significance at α = 5%."),
          tags$li(paste0("Trait '", best, "' has the smallest CD (", round(min(rv$cd$CD_5pct),3), ") — most sensitive to treatment differences.")),
          tags$li("Smaller SEd → better experimental precision and more discriminating power."),
          tags$li("Compare observed mean differences with CD: if diff > CD, treatments are significantly different.")
        )
    )
  })
  
  # ── G×E ANOVA ────────────────────────────────────────────────────────────────
  observeEvent(input$run_gxe, {
    req(rv$df, input$gxe_trait)
    rv$gxe <- run_gxe_anova(rv$df, input$gxe_trait)
    rv$selected_gxe_trait <- input$gxe_trait
  })
  
  output$gxe_table <- renderDT({
    req(rv$gxe)
    datatable(rv$gxe, options=list(pageLength=10), rownames=FALSE) %>%
      formatStyle("Signif", color=styleEqual(c("***","**","*"),c("#f85149","#ff9f1c","#3fb950")))
  })
  
  gxe_plot_obj <- reactive({
    req(rv$df, rv$selected_gxe_trait)
    tr <- rv$selected_gxe_trait
    means <- rv$df %>%
      group_by(GEN, ENV) %>%
      summarise(val=mean(.data[[tr]], na.rm=TRUE), .groups="drop")
    ggplot(means, aes(x=ENV, y=val, color=GEN, group=GEN)) +
      geom_line(linewidth=1.1, alpha=0.9) +
      geom_point(size=3) +
      labs(title=paste("G×E Interaction:", tr), x="Environment", y=paste("Mean", tr)) +
      theme_minimal(base_size=13) +
      theme(
        plot.background=element_rect(fill="#161b22",color=NA),
        panel.background=element_rect(fill="#161b22",color=NA),
        panel.grid.major=element_line(color="#21262d"),
        panel.grid.minor=element_blank(),
        axis.text=element_text(color="#8b949e"),
        axis.title=element_text(color="#8b949e"),
        plot.title=element_text(color="#e6edf3",face="bold",size=14),
        legend.background=element_rect(fill="#161b22",color=NA),
        legend.text=element_text(color="#8b949e"),
        legend.title=element_text(color="#8b949e")
      )
  })
  
  output$gxe_plot <- renderPlot({ gxe_plot_obj() })
  
  output$gxe_interp <- renderUI({
    req(rv$gxe)
    gxe_row <- rv$gxe[grepl("GEN:ENV|GEN.ENV", rv$gxe$Source), ]
    p_gxe   <- if(nrow(gxe_row)>0) gxe_row$P_value[1] else NA
    div(class="interp-box",
        tags$h4("💡 Interpretation"),
        tags$ul(
          tags$li(paste0("G×E Interaction p-value: ", if(is.na(p_gxe)) "N/A" else p_gxe,
                         if(!is.na(p_gxe) && p_gxe < 0.05) " *** Significant" else " (not significant)")),
          tags$li(if(!is.na(p_gxe) && p_gxe < 0.05)
            "⚠️ Significant G×E: genotype rankings are NOT consistent across environments. Stability analysis is strongly recommended."
            else
              "✅ Non-significant G×E: genotype rankings are relatively consistent — a single best genotype may exist."),
          tags$li("Crossed lines in the interaction plot visually indicate rank changes (crossover interactions).")
        )
    )
  })
  
  # ── ENV ANOVA ────────────────────────────────────────────────────────────────
  env_res_list <- reactiveVal(NULL)
  
  observeEvent(input$run_env, {
    req(rv$df, rv$traits)
    res <- lapply(rv$traits, function(tr) run_env_anova(rv$df, tr))
    names(res) <- rv$traits
    env_res_list(res)
  })
  
  output$env_trait_tabs <- renderUI({
    req(env_res_list())
    tabs <- lapply(names(env_res_list()), function(tr) {
      tabPanel(tr, DTOutput(paste0("env_tbl_", tr)))
    })
    do.call(tabsetPanel, tabs)
  })
  
  observe({
    req(env_res_list())
    lapply(names(env_res_list()), function(tr) {
      local({
        tr_local <- tr
        output[[paste0("env_tbl_", tr_local)]] <- renderDT({
          datatable(env_res_list()[[tr_local]], options=list(pageLength=5), rownames=FALSE) %>%
            formatStyle("Signif", color=styleEqual(c("***","**","*"),c("#f85149","#ff9f1c","#3fb950")))
        })
      })
    })
  })
  
  output$env_interp <- renderUI({
    req(env_res_list())
    div(class="interp-box",
        tags$h4("💡 Interpretation"),
        tags$ul(
          tags$li("Significant ENV effect (p < 0.05) confirms that environments differ meaningfully in productivity."),
          tags$li("This justifies conducting multi-environment trials and computing stability statistics."),
          tags$li("Non-significant ENV effect would suggest all environments provide similar conditions.")
        )
    )
  })
  
  # ── STABILITY ────────────────────────────────────────────────────────────────
  observeEvent(input$run_stab, {
    req(rv$df, input$stab_trait)
    out <- eberhart_stability(rv$df, input$stab_trait)
    rv$stab          <- out$stability
    rv$stab_env      <- out$env_idx
    rv$stab_reg      <- out$data_reg
    rv$selected_stab_trait <- input$stab_trait
  })
  
  output$stab_table <- renderDT({
    req(rv$stab)
    datatable(rv$stab, options=list(pageLength=15), rownames=FALSE) %>%
      formatRound(c("Mean","bi","s2d"), 4) %>%
      formatStyle("Stability_Class",
                  color = styleEqual(
                    c("✅ Stable","📈 Responsive","📉 Below avg env","⚠️ Unstable"),
                    c("#3fb950","#58a6ff","#ff9f1c","#f85149")
                  ))
  })
  
  stab_plot_obj <- reactive({
    req(rv$stab)
    grand <- mean(rv$stab$Mean)
    ggplot(rv$stab, aes(x=Mean, y=bi, color=Stability_Class, label=GEN)) +
      geom_hline(yintercept=1,  linetype="dashed", color="#8b949e", linewidth=0.8) +
      geom_vline(xintercept=grand, linetype="dashed", color="#8b949e", linewidth=0.8) +
      geom_point(size=5, alpha=0.9) +
      geom_text(vjust=-0.9, size=3.8, fontface="bold") +
      scale_color_manual(values=c("✅ Stable"="#3fb950","📈 Responsive"="#58a6ff","📉 Below avg env"="#ff9f1c","⚠️ Unstable"="#f85149")) +
      labs(title=paste("Stability Biplot:", rv$selected_stab_trait),
           subtitle="Dashed lines: grand mean (vertical) and bi = 1 (horizontal)",
           x="Mean Yield", y="Regression Coefficient (bi)", color="Class") +
      theme_minimal(base_size=13) +
      theme(
        plot.background=element_rect(fill="#161b22",color=NA),
        panel.background=element_rect(fill="#161b22",color=NA),
        panel.grid.major=element_line(color="#21262d"),
        panel.grid.minor=element_blank(),
        axis.text=element_text(color="#8b949e"),
        axis.title=element_text(color="#8b949e"),
        plot.title=element_text(color="#e6edf3",face="bold",size=14),
        plot.subtitle=element_text(color="#8b949e",size=11),
        legend.background=element_rect(fill="#161b22",color=NA),
        legend.text=element_text(color="#8b949e",size=10),
        legend.title=element_text(color="#8b949e",size=10)
      )
  })
  
  output$stab_plot <- renderPlot({ stab_plot_obj() })
  
  stab_reg_obj <- reactive({
    req(rv$stab_reg, rv$stab_env)
    dr <- rv$stab_reg %>% left_join(rv$stab_env, by="ENV")
    ggplot(dr, aes(x=Ij, y=mean_val, color=GEN, group=GEN)) +
      geom_point(size=2.5, alpha=0.7) +
      geom_smooth(method="lm", se=FALSE, linewidth=1) +
      labs(title=paste("Regression of Genotype Means on Env. Index:", rv$selected_stab_trait),
           x="Environment Index (Ij)", y="Genotype Mean") +
      theme_minimal(base_size=13) +
      theme(
        plot.background=element_rect(fill="#161b22",color=NA),
        panel.background=element_rect(fill="#161b22",color=NA),
        panel.grid.major=element_line(color="#21262d"),
        panel.grid.minor=element_blank(),
        axis.text=element_text(color="#8b949e"),
        axis.title=element_text(color="#8b949e"),
        plot.title=element_text(color="#e6edf3",face="bold",size=13),
        legend.background=element_rect(fill="#161b22",color=NA),
        legend.text=element_text(color="#8b949e")
      )
  })
  
  output$stab_regplot <- renderPlot({ stab_reg_obj() })
  
  output$stab_interp <- renderUI({
    req(rv$stab)
    best_row <- rv$stab %>% filter(grepl("Stable", Stability_Class)) %>% arrange(desc(Mean))
    div(class="interp-box",
        tags$h4("💡 Interpretation"),
        tags$ul(
          tags$li(if(nrow(best_row)>0) paste0("Best stable genotype: ", best_row$GEN[1], " (Mean = ", best_row$Mean[1], ", bi = ", best_row$bi[1], ", s²d = ", best_row$s2d[1], ")") else "No genotype classified as stable. Consider selecting high-mean or specific adaptation genotypes."),
          tags$li("Genotypes in upper-left quadrant (high mean, bi ≈ 1): ideal for broad adaptation."),
          tags$li("Genotypes in upper-right (high mean, bi > 1): suited for favorable environments."),
          tags$li("Steeper regression lines in regression plot = more responsive genotypes.")
        )
    )
  })
  
  # ── Downloads ─────────────────────────────────────────────────────────────────
  make_dl_csv <- function(data_rv) {
    downloadHandler(
      filename = function() paste0(deparse(substitute(data_rv)), "_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(data_rv())
        write.csv(data_rv(), file, row.names=FALSE)
      }
    )
  }
  
  output$dl_data_csv <- downloadHandler(
    filename = function() paste0("MET_data_", Sys.Date(), ".csv"),
    content  = function(file) { req(rv$df); write.csv(rv$df, file, row.names=FALSE) }
  )
  
  output$dl_cd_csv <- downloadHandler(
    filename = function() paste0("CD_results_", Sys.Date(), ".csv"),
    content  = function(file) { req(rv$cd); write.csv(rv$cd, file, row.names=FALSE) }
  )
  
  output$dl_gxe_csv <- downloadHandler(
    filename = function() paste0("GxE_ANOVA_", Sys.Date(), ".csv"),
    content  = function(file) { req(rv$gxe); write.csv(rv$gxe, file, row.names=FALSE) }
  )
  
  output$dl_env_csv <- downloadHandler(
    filename = function() paste0("Env_ANOVA_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(env_res_list())
      out <- bind_rows(lapply(names(env_res_list()), function(tr) {
        cbind(Trait=tr, env_res_list()[[tr]])
      }))
      write.csv(out, file, row.names=FALSE)
    }
  )
  
  output$dl_stab_csv <- downloadHandler(
    filename = function() paste0("Stability_", Sys.Date(), ".csv"),
    content  = function(file) { req(rv$stab); write.csv(rv$stab, file, row.names=FALSE) }
  )
  
  output$dl_all_txt <- downloadHandler(
    filename = function() paste0("All_Results_", Sys.Date(), ".txt"),
    content  = function(file) {
      sink(file)
      cat("========== CRITICAL DIFFERENCE (CD) ==========\n\n")
      if(!is.null(rv$cd)) print(rv$cd, row.names=FALSE)
      cat("\n\n========== G × E ANOVA ==========\n")
      if(!is.null(rv$gxe)) { cat("\n--- Trait:", rv$selected_gxe_trait, "---\n"); print(rv$gxe, row.names=FALSE) }
      cat("\n\n========== ENVIRONMENT-ONLY ANOVA ==========\n")
      if(!is.null(env_res_list())) {
        lapply(names(env_res_list()), function(tr) {
          cat("\n--- Trait:", tr, "---\n"); print(env_res_list()[[tr]], row.names=FALSE)
        })
      }
      cat("\n\n========== EBERHART-RUSSELL STABILITY ==========\n")
      if(!is.null(rv$stab)) { cat("\n--- Trait:", rv$selected_stab_trait, "---\n"); print(rv$stab, row.names=FALSE) }
      sink()
    }
  )
  
  plot_dl <- function(plot_fn, fname) {
    downloadHandler(
      filename = function() paste0(fname, "_", Sys.Date(), ".png"),
      content  = function(file) {
        p <- plot_fn()
        ggsave(file, plot=p, width=10, height=6, dpi=200, bg="#161b22")
      }
    )
  }
  
  output$dl_cd_plot   <- plot_dl(cd_plot_obj,   "CD_plot")
  output$dl_stab_plot <- plot_dl(stab_plot_obj, "Stability_biplot")
  output$dl_reg_plot  <- plot_dl(stab_reg_obj,  "Regression_plot")
  output$dl_gxe_plot  <- plot_dl(gxe_plot_obj,  "GxE_interaction_plot")
}

shinyApp(ui, server)

