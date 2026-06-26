###############################################################################
# 38_si_nonlinearity.R
#
# Standalone non-linearity check for the four headline ENSO IRFs (extracted
# from analysis notebooks/Extended Data.Rmd). Each panel compares linear vs
# polynomial (quadratic, cubic) NINO3.4 specifications.
#
# Panel C (WR teleconnected yields) uses the v2 panel + gold-ticket #3 spec:
#   outcome = log(yield/sd(yield)) (% harvest)
#   FE      = VarLocationGrain[Year]
#   controls = i(decade) + l(nino34, 1:3)
# All other panels reproduce the original notebook spec.
#
# Output:
#   analysis/output/figures/appendix/FigA_nonlin.pdf
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
})
setFixest_notes(FALSE)

.script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg)))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                error = function(e) getwd())
}
SCRIPT_DIR <- .script_dir()
ROOT       <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
DATA_PROC  <- file.path(ROOT, "processed data")
fig_out    <- file.path(SCRIPT_DIR, "output", "figures", "appendix")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

# ── Data ─────────────────────────────────────────────────────────────────────
famines <- read.csv(file.path(DATA_PROC, "famine_region_data.csv")) |>
  mutate(Decade = floor(Year / 10) * 10)

price <- read.csv(file.path(DATA_PROC, "price_2023_enso.csv"))
if (!"Decade" %in% names(price)) price$Decade <- floor(price$Year / 10) * 10

fishprice <- read.csv(file.path(DATA_PROC, "fishprice_enso.csv"))
if (!"Decade" %in% names(fishprice)) fishprice$Decade <- floor(fishprice$Year / 10) * 10

# v2 yield panel, raw logyield outcome (β reads as % harvest).
yield_v2 <- read.csv(file.path(DATA_PROC, "yield_ljungqvist_v2.csv")) |>
  mutate(decade = floor(Year / 10) * 10,
         Year2  = Year^2,
         Year3  = Year^3) |>
  filter(!is.na(yield), yield > 0)

# ── Helper: predict response on grid of nino34 values ───────────────────────
predict_poly <- function(model, var = "nino34",
                          xmin = -1, xmax = 1.5, ref = 0,
                          step = 0.01) {
  grid <- data.frame(x = seq(xmin, xmax, by = step))
  cn   <- names(coef(model))
  V    <- vcov(model)
  # locate polynomial terms in the form poly(var, k, raw=TRUE)k OR plain var
  poly_terms <- grep(paste0("poly\\(", var), cn, value = TRUE)
  plain_term <- if (var %in% cn) var else NULL
  if (length(poly_terms) > 0) {
    # extract degree from the term name: ...)1, ...)2, ...)3
    degs <- as.integer(sub(".*\\)([0-9]+)$", "\\1", poly_terms))
    coefs <- coef(model)[poly_terms]
    # response at x (relative to reference)
    f_x <- function(x) sum(coefs * (x^degs - ref^degs))
    g_x <- function(x) (x^degs - ref^degs)  # gradient vector
    grid$response <- vapply(grid$x, f_x, numeric(1))
    se <- vapply(grid$x, function(x) {
      g <- g_x(x); sqrt(as.numeric(t(g) %*% V[poly_terms, poly_terms] %*% g))
    }, numeric(1))
  } else if (!is.null(plain_term)) {
    b <- coef(model)[plain_term]
    s <- sqrt(V[plain_term, plain_term])
    grid$response <- b * (grid$x - ref)
    se <- abs(grid$x - ref) * s
  } else {
    stop("Could not find ", var, " in model coefficients")
  }
  grid$lower <- grid$response - 1.96 * se
  grid$upper <- grid$response + 1.96 * se
  setNames(grid, c("nino34", "response", "lower", "upper"))
}

pp <- function(m) predict_poly(m, "nino34", xmin = -1, xmax = 1.5, ref = 0)

make_panel <- function(lin, p2, p3, y_lab, tag) {
  spec_col <- c("Linear" = "#4477AA",
                "Quadratic" = "tomato",
                "Cubic" = "orange")
  df_lin <- lin |> mutate(spec = "Linear")
  df_p2  <- p2  |> mutate(spec = "Quadratic")
  df_p3  <- p3  |> mutate(spec = "Cubic")
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(data = df_lin,
                aes(x = nino34, ymin = lower, ymax = upper),
                fill = spec_col["Linear"], alpha = 0.2) +
    geom_line(data = bind_rows(df_lin, df_p2, df_p3),
              aes(x = nino34, y = response, colour = spec),
              linewidth = 0.9) +
    scale_colour_manual(values = spec_col, name = NULL,
                        guide = guide_legend(
                          override.aes = list(linewidth = 1.5))) +
    labs(x = "NINO3.4", y = y_lab, tag = tag) +
    theme_classic(base_size = 13) +
    theme(legend.position  = "none",
          plot.tag         = element_text(size = 18, face = "bold"),
          axis.title.y     = element_text(size = 10),
          panel.grid.minor = element_blank())
}

# ── Panel A: famine onset, Central Europe ──────────────────────────────────
ce <- famines |> filter(Region == "Central Europe")
m1_lin <- feols(Famine_start ~ nino34 | Region,
                data = ce, panel.id = c("Region","Year"), vcov = DK ~ Year)
m1_p2  <- feols(Famine_start ~ poly(nino34, 2, raw = TRUE) | Region,
                data = ce, panel.id = c("Region","Year"), vcov = DK ~ Year)
m1_p3  <- feols(Famine_start ~ poly(nino34, 3, raw = TRUE) | Region,
                data = ce, panel.id = c("Region","Year"), vcov = DK ~ Year)
p1 <- make_panel(pp(m1_lin), pp(m1_p2), pp(m1_p3),
                 y_lab = "P(famine onset)", tag = "a")

# ── Panel B: grain prices, h = 1 ────────────────────────────────────────────
m2_lin <- feols(f(logprice, 1) - l(logprice, 1) ~ nino34 + i(Decade) | Location,
                data = price, panel.id = c("Location","Year"), vcov = DK ~ Year)
m2_p2  <- feols(f(logprice, 1) - l(logprice, 1) ~ poly(nino34, 2, raw = TRUE) +
                  i(Decade) | Location,
                data = price, panel.id = c("Location","Year"), vcov = DK ~ Year)
m2_p3  <- feols(f(logprice, 1) - l(logprice, 1) ~ poly(nino34, 3, raw = TRUE) +
                  i(Decade) | Location,
                data = price, panel.id = c("Location","Year"), vcov = DK ~ Year)
p2 <- make_panel(pp(m2_lin), pp(m2_p2), pp(m2_p3),
                 y_lab = "Log grain price change (h = 1)", tag = "b")

# ── Panel C: WR teleconnected yields, h = 0 (v2 + GT#3) ─────────────────────
WR_tele <- yield_v2 |> filter(Grain %in% c("Wheat","Rye"),
                               teleco_PDSI_10 == 1)
m3_lin <- feols(f(logyield, 0) - l(logyield, 1) ~ nino34 + l(nino34, 1:3) +
                  i(decade) | VarLocationGrain[Year],
                data = WR_tele, panel.id = c("VarLocationGrain","Year"),
                vcov = DK ~ Year)
m3_p2  <- feols(f(logyield, 0) - l(logyield, 1) ~ poly(nino34, 2, raw = TRUE) +
                  l(nino34, 1:3) + i(decade) |
                  VarLocationGrain[Year],
                data = WR_tele, panel.id = c("VarLocationGrain","Year"),
                vcov = DK ~ Year)
m3_p3  <- feols(f(logyield, 0) - l(logyield, 1) ~ poly(nino34, 3, raw = TRUE) +
                  l(nino34, 1:3) + i(decade) |
                  VarLocationGrain[Year],
                data = WR_tele, panel.id = c("VarLocationGrain","Year"),
                vcov = DK ~ Year)
p3 <- make_panel(pp(m3_lin), pp(m3_p2), pp(m3_p3),
                 y_lab = "Log yield change (h = 0)", tag = "c")

# ── Panel D: fish prices, h = 3 ─────────────────────────────────────────────
m4_lin <- feols(f(logprice, 3) - l(logprice, 1) ~ nino34 + i(Decade) | LocationSpecies,
                data = fishprice, panel.id = c("LocationSpecies","Year"),
                vcov = DK ~ Year)
m4_p2  <- feols(f(logprice, 3) - l(logprice, 1) ~ poly(nino34, 2, raw = TRUE) +
                  i(Decade) | LocationSpecies,
                data = fishprice, panel.id = c("LocationSpecies","Year"),
                vcov = DK ~ Year)
m4_p3  <- feols(f(logprice, 3) - l(logprice, 1) ~ poly(nino34, 3, raw = TRUE) +
                  i(Decade) | LocationSpecies,
                data = fishprice, panel.id = c("LocationSpecies","Year"),
                vcov = DK ~ Year)
p4 <- make_panel(pp(m4_lin), pp(m4_p2), pp(m4_p3),
                 y_lab = "Log fish price change (h = 3)", tag = "d")

# ── Assemble ────────────────────────────────────────────────────────────────
p_legend <- p1 +
  theme(legend.position  = "top",
        legend.text      = element_text(size = 14),
        legend.key.width = unit(1.5, "cm"),
        legend.spacing.x = unit(0.4, "cm"))
legend_grob <- ggpubr::get_legend(p_legend)
p_grid <- (p1 | p2) / (p3 | p4)
p_all  <- ggpubr::ggarrange(
  ggpubr::as_ggplot(legend_grob),
  p_grid,
  ncol    = 1,
  heights = c(0.06, 1)
)

out_path <- file.path(fig_out, "FigA_nonlin.pdf")
ggsave(out_path, p_all, width = 8, height = 6, device = cairo_pdf)
message("Saved ", out_path)
