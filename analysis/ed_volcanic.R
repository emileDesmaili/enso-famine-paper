# Auto-converted from analysis notebooks/volcanic_forcing.Rmd (do not edit the Rmd —
# edit this script or re-purl if the Rmd is kept in the archive).
# Produces: analysis/output/figures/extended data/figED_volcanic.pdf

.script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg)))
  } else {
    tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
             error = function(e) getwd())
  }
}
SCRIPT_DIR <- .script_dir()
ROOT_DIR   <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(fixest)
  library(readxl)
  library(purrr)
  library(broom)
  library(scales)
})

source(file.path(ROOT_DIR, "emileRegs.R"))
setFixest_notes(FALSE)

# ── Colour palette ───────────────────────────────────────────────────────────
# Panel (a): vibrant, distinct hues for the four time series.
COL_ENSO  <- "#1F4E9E"          # vivid royal blue  – ENSO series
COL_GVF   <- "#C0392B"          # deep crimson      – GVF series
COL_GR    <- "#0F8B8D"          # teal              – Greenland sulfate
COL_ANT   <- "#F1A208"          # warm gold         – Antarctica sulfate
# Panel (b): black baseline + coral V.F.-controlled overlay.
COL_BASE  <- "black"            # black             – baseline price IRF
COL_VOLC  <- "#D62828"          # vibrant firebrick – V.F.-controlled IRF
# Panel (c): deep blue baseline + lighter blue V.F.-controlled.
COL_C_BASE <- "#6495ED"         # cornflower blue   – panel-c baseline
COL_C_VOLC <- "#D62828"         # vibrant firebrick – panel-c V.F.-controlled

nature_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line         = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length = unit(3, "pt"),
      axis.text         = element_text(size = base_size - 1, colour = "black",
                                       face = "bold"),
      axis.title        = element_text(size = base_size, colour = "black",
                                       face = "bold"),
      legend.position   = "bottom",
      legend.key.size   = unit(10, "pt"),
      legend.text       = element_text(size = base_size - 1, face = "bold"),
      legend.title      = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      legend.spacing.x  = unit(5, "pt"),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = base_size),
      plot.tag          = element_text(face = "bold", size = base_size + 4,
                                       colour = "black"),
      plot.tag.position = "topleft",
      plot.margin       = margin(6, 8, 6, 8, "pt"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
}

raw <- read_excel(file.path(ROOT_DIR, "data/volcanic_data.xlsx"),
                  sheet = "1 - Volcanic_Forcing",
                  col_names = FALSE)

# Row 1 = header level 1, row 2 = header level 2, rows 3+ = data
# Columns (1-indexed in R):
#   1  Year
#   2  Source (1=tropical, 2=NH, 3=SH)
#   3  Greenland sulfate [kg km-2]
#   4  Greenland uncertainty
#   5  Antarctica sulfate [kg km-2]
#   6  Antarctica uncertainty
#   7  Global Volcanic Forcing GVF [W m-2]
#   8  GVF uncertainty

volc_raw <- raw[-c(1, 2), 1:8]
colnames(volc_raw) <- c(
  "Year", "Source",
  "greenland_sulfate", "greenland_unc",
  "antarctica_sulfate", "antarctica_unc",
  "GVF", "GVF_unc"
)

volc_raw <- volc_raw %>%
  mutate(across(everything(), as.numeric)) %>%
  filter(!is.na(Year), Year >= 1500, Year <= 1900) %>%
  arrange(Year)

cat("Volcanic events in sample (1500-1900):", nrow(volc_raw), "\n")
print(volc_raw)

year_grid <- tibble(Year = 1500:1900)

volc <- year_grid %>%
  left_join(volc_raw, by = "Year") %>%
  mutate(
    GVF                = replace_na(GVF, 0),
    greenland_sulfate  = replace_na(greenland_sulfate, 0),
    antarctica_sulfate = replace_na(antarctica_sulfate, 0),
    eruption           = as.integer(GVF != 0)
  ) %>%
  dplyr::select(Year, GVF, greenland_sulfate, antarctica_sulfate, eruption) %>%
  mutate(
    GVF_l1 = lag(GVF, 1),
    GVF_l2 = lag(GVF, 2),
    GVF_l3 = lag(GVF, 3)
  )

summary(volc[, c("GVF", "greenland_sulfate", "antarctica_sulfate")])

price <- read.csv(file.path(ROOT_DIR, "processed data/price_2023_enso.csv")) %>%
  mutate(Decade = floor(Year / 10) * 10) %>%
  drop_na() %>% distinct() %>% arrange(Location, Year)

fishprice <- read.csv(file.path(ROOT_DIR, "processed data/fishprice_enso.csv")) %>%
  mutate(Decade = floor(Year / 10) * 10) %>%
  arrange(LocationSpecies, Year)

famines <- read.csv(file.path(ROOT_DIR, "processed data/famine_region_data.csv")) %>%
  mutate(Decade = floor(Year / 10) * 10) %>%
  arrange(Region, Year)

price     <- left_join(price,     volc, by = "Year")
fishprice <- left_join(fishprice, volc, by = "Year")
famines   <- left_join(famines,   volc, by = "Year")

cat("Grain price obs after merge:", nrow(price), "\n")
cat("Fish price obs after merge: ", nrow(fishprice), "\n")
cat("Famine obs after merge:     ", nrow(famines), "\n")

ts <- price %>%
  group_by(Year) %>%
  summarise(
    ENSO               = mean(nino34,             na.rm = TRUE),
    GVF                = mean(GVF,                na.rm = TRUE),
    greenland_sulfate  = mean(greenland_sulfate,  na.rm = TRUE),
    antarctica_sulfate = mean(antarctica_sulfate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(ENSO, GVF, greenland_sulfate, antarctica_sulfate),
                ~ as.numeric(scale(.x))))

fmt_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

cor_gvf  <- cor.test(ts$ENSO, ts$GVF)
cor_gr   <- cor.test(ts$ENSO, ts$greenland_sulfate)
cor_ant  <- cor.test(ts$ENSO, ts$antarctica_sulfate)

cat(sprintf("ENSO vs GVF:                r = %.3f, %s\n",
            cor_gvf$estimate,  fmt_p(cor_gvf$p.value)))
cat(sprintf("ENSO vs Greenland sulfate:  r = %.3f, %s\n",
            cor_gr$estimate,   fmt_p(cor_gr$p.value)))
cat(sprintf("ENSO vs Antarctica sulfate: r = %.3f, %s\n",
            cor_ant$estimate,  fmt_p(cor_ant$p.value)))

run_lp <- function(df, controls, fe, panel_id) {
  lp_panel(
    data         = df,
    outcome      = "logprice",
    main_var     = "nino34",
    controls     = controls,
    horizon      = 10,
    fe           = fe,
    panel_id     = panel_id,
    vcov_formula = DK ~ Year
  )
}

irf_grain_base <- run_lp(
  price,
  controls = "i(Decade) + l(nino34, 1:3)",
  fe       = "Location",
  panel_id = c("Location", "Year")
) %>% mutate(spec = "Baseline")

irf_grain_volc <- run_lp(
  price %>% drop_na(GVF),
  controls = "i(Decade) + GVF + greenland_sulfate + antarctica_sulfate + l(nino34, 1:3)",
  fe       = "Location",
  panel_id = c("Location", "Year")
) %>% mutate(spec = "Controlling for V.F.")

irf_grain <- bind_rows(irf_grain_base, irf_grain_volc) %>%
  mutate(spec = factor(spec, levels = c("Baseline", "Controlling for V.F.")))

m_base <- feols(
  Famine_start ~ nino34:i(Region) | Region + Decade,
  data     = famines,
  panel.id = c("Region", "Year"),
  vcov     = DK ~ Year
)

m_volc <- feols(
  Famine_start ~ nino34:i(Region) + GVF + greenland_sulfate + antarctica_sulfate | Region + Decade,
  data     = famines %>% drop_na(GVF),
  panel.id = c("Region", "Year"),
  vcov     = DK ~ Year
)

etable(m_base, m_volc,
       dict    = c("nino34" = "NINO3.4", "GVF" = "GVF T",
                   "GVF_l1" = "GVF T-1", "GVF_l2" = "GVF T-2",
                   "GVF_l3" = "GVF T-3"),
       fitstat = c("n", "r2"))

# Extract region coefficients with 90 and 95 % CIs
extract_region_coefs <- function(mod, spec_label) {
  b   <- coef(mod)
  V   <- vcov(mod)
  nms <- names(b)
  idx <- grep("nino34", nms)

  tibble(
    term     = nms[idx],
    estimate = b[idx],
    se       = sqrt(diag(V)[idx]),
    spec     = spec_label
  ) %>%
    mutate(
      lo95 = estimate - 1.960 * se,
      hi95 = estimate + 1.960 * se,
      lo90 = estimate - 1.645 * se,
      hi90 = estimate + 1.645 * se,
      Region = term %>%
        gsub("nino34:Region::", "", .) %>%
        gsub("nino34:i\\(Region\\)::", "", .) %>%
        trimws()
    )
}

coef_base <- extract_region_coefs(m_base, "Baseline")
coef_volc <- extract_region_coefs(m_volc, "Controlling for V.F.")
coef_df   <- bind_rows(coef_base, coef_volc) %>%
  mutate(spec = factor(spec, levels = c("Baseline", "Controlling for V.F.")))

print(coef_df[, c("Region", "spec", "estimate", "lo95", "hi95")])

# Build annotation labels
make_cor_label <- function(cor_obj, name) {
  sprintf("%s:  r = %.2f  (%s)", name,
          cor_obj$estimate, fmt_p(cor_obj$p.value))
}

ann_df <- tibble(
  label  = c("ENSO (NINO3.4)",
             make_cor_label(cor_gvf, "GVF"),
             make_cor_label(cor_gr,  "Greenland"),
             make_cor_label(cor_ant, "Antarctica")),
  colour = c(COL_ENSO, COL_GVF, COL_GR, COL_ANT)
)

ts_long <- ts %>%
  pivot_longer(-Year, names_to = "Index", values_to = "value") %>%
  mutate(Index = factor(Index,
    levels = c("ENSO", "GVF", "greenland_sulfate", "antarctica_sulfate"),
    labels = c("ENSO (NINO3.4)", "Global V.F. (GVF)",
               "Greenland sulfate", "Antarctica sulfate")))

panel_a <- ggplot(ts_long, aes(Year, value, colour = Index, linewidth = Index)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3, linetype = "dashed") +
  geom_line(alpha = 0.88) +
  # Background label box for readability
  annotate("rect",
           xmin = 1500, xmax = 1575,
           ymin = max(ts_long$value, na.rm = TRUE) * 0.97 -
                  diff(range(ts_long$value, na.rm = TRUE)) * 0.28,
           ymax = max(ts_long$value, na.rm = TRUE) * 0.97,
           fill = "white", colour = NA) +
  # One annotation per line, coloured by series
  annotate("text",
           x = 1505,
           y = max(ts_long$value, na.rm = TRUE) * 0.97 -
               diff(range(ts_long$value, na.rm = TRUE)) *
               c(0.04, 0.10, 0.16, 0.22),
           label    = ann_df$label,
           colour   = ann_df$colour,
           hjust    = 0, vjust = 0.5,
           size     = 4, fontface = "bold") +
  scale_colour_manual(
    values = c("ENSO (NINO3.4)"     = COL_ENSO,
               "Global V.F. (GVF)"  = COL_GVF,
               "Greenland sulfate"  = COL_GR,
               "Antarctica sulfate" = COL_ANT)
  ) +
  scale_linewidth_manual(
    values = c("ENSO (NINO3.4)"     = 0.9,
               "Global V.F. (GVF)"  = 0.7,
               "Greenland sulfate"  = 0.7,
               "Antarctica sulfate" = 0.7)
  ) +
  guides(linewidth = "none", colour = "none") +
  scale_x_continuous(breaks = seq(1500, 1900, 50)) +
  labs(x = "Year", y = "Standardised value",
       title = "ENSO and volcanic forcing") +
  nature_theme() +
  theme(legend.position = "none")

panel_a

# Baseline: solid line with CI ribbon; V.F.-controlled: dashed overlay, no CI.
irf_base_only <- dplyr::filter(irf_grain, spec == "Baseline")
irf_volc_only <- dplyr::filter(irf_grain, spec == "Controlling for V.F.")

panel_b <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 1.0,
             linetype = "dashed") +
  # baseline 95% CI ribbon (matches fig 5 style: cornflowerblue, alpha 0.20)
  geom_ribbon(data = irf_base_only,
              aes(x = horizon, ymin = irf_down, ymax = irf_up),
              fill = "cornflowerblue", alpha = 0.20, colour = NA) +
  # baseline line
  geom_line(data = irf_base_only,
            aes(x = horizon, y = irf_mean, colour = spec, linetype = spec),
            linewidth = 1.6) +
  # controls-for-V.F. line only
  geom_line(data = irf_volc_only,
            aes(x = horizon, y = irf_mean, colour = spec, linetype = spec),
            linewidth = 1.6) +
  scale_x_continuous(breaks = seq(0, 10, 2), expand = c(0.02, 0)) +
  scale_colour_manual(values   = c("Baseline" = COL_BASE,
                                    "Controlling for V.F." = COL_VOLC)) +
  scale_linetype_manual(values = c("Baseline" = "solid",
                                    "Controlling for V.F." = "longdash")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  guides(colour   = guide_legend(override.aes = list(linewidth = 1.6)),
         linetype = guide_legend(override.aes = list(linewidth = 1.6))) +
  labs(x = "Horizon (years)", y = "% Grain price response",
       title = "ENSO effect on grain prices") +
  nature_theme()

panel_b

# Dodge position for two specs
pd <- position_dodge(width = 0.55)

panel_c <- ggplot(coef_df,
                  aes(x = Region, y = estimate,
                      colour = spec, shape = spec)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35,
             linetype = "dashed") +
  # 95 % CI (thin whisker)
  geom_errorbar(aes(ymin = lo95, ymax = hi95),
                width = 0.0, linewidth = 0.55,
                position = pd) +
  # 90 % CI (thick whisker)
  geom_errorbar(aes(ymin = lo90, ymax = hi90),
                width = 0.0, linewidth = 1.5,
                position = pd) +
  geom_point(size = 2.8, position = pd) +
  # Highlight Central Europe with a light grey-blue shaded band
  annotate("rect",
           xmin = 0.5, xmax = 1.5,
           ymin = -Inf, ymax = Inf,
           fill = "#B7C9D8", alpha = 0.25) +
  scale_colour_manual(values = c(COL_C_BASE, COL_C_VOLC)) +
  scale_shape_manual(values  = c(16, 17)) +
  guides(colour = guide_legend(override.aes = list(size = 2.5)),
         shape  = guide_legend(override.aes = list(size = 2.5))) +
  labs(x = NULL, y = "ENSO effect on famine onset prob.",
       title = "Famine onset by region") +
  nature_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10))

panel_c

fig <- panel_a / ((panel_b | panel_c) + plot_layout(widths = c(0.85, 1.25))) +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag   = element_text(face = "bold", size = 20, colour = "black"),
        plot.title = element_text(face = "bold", size = 15))

fig

# IRF of GVF on grain prices controlling for ENSO
irf_gvf_grain <- lp_panel(
  data         = price,
  outcome      = "logprice",
  main_var     = "GVF",
  controls     = "i(Decade) + nino34",
  horizon      = 10,
  fe           = "Location",
  panel_id     = c("Location", "Year"),
  vcov_formula = DK ~ Year
) %>% mutate(forcing = "Global V.F. (GVF)")

# Also run for Greenland and Antarctica separately
irf_gvf_gr <- lp_panel(
  data         = price,
  outcome      = "logprice",
  main_var     = "greenland_sulfate",
  controls     = "i(Decade) + nino34",
  horizon      = 10,
  fe           = "Location",
  panel_id     = c("Location", "Year"),
  vcov_formula = DK ~ Year
) %>% mutate(forcing = "Greenland sulfate")

irf_gvf_ant <- lp_panel(
  data         = price,
  outcome      = "logprice",
  main_var     = "antarctica_sulfate",
  controls     = "i(Decade) + nino34",
  horizon      = 10,
  fe           = "Location",
  panel_id     = c("Location", "Year"),
  vcov_formula = DK ~ Year
) %>% mutate(forcing = "Antarctica sulfate")

irf_volc_all <- bind_rows(irf_gvf_grain, irf_gvf_gr, irf_gvf_ant) %>%
  mutate(forcing = factor(forcing,
    levels = c("Global V.F. (GVF)", "Greenland sulfate", "Antarctica sulfate")))

COL_VOLC3 <- c("Global V.F. (GVF)"  = COL_GVF,
               "Greenland sulfate"  = COL_GR,
               "Antarctica sulfate" = COL_ANT)

panel_d <- ggplot(irf_volc_all,
                  aes(horizon, irf_mean,
                      colour = forcing, fill = forcing)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.0) +
  scale_x_continuous(breaks = seq(0, 10, 2), expand = c(0.02, 0)) +
  scale_colour_manual(values = COL_VOLC3) +
  scale_fill_manual(values   = COL_VOLC3) +
  scale_y_continuous(labels  = percent_format(accuracy = 0.1)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.4)),
         fill   = guide_legend(override.aes = list(alpha = 0.2))) +
  labs(x = "Horizon (years)", y = "% Grain price response",
       title = "V.F. effect on grain prices (controlling for ENSO)") +
  nature_theme()

panel_d

m_volc_region <- feols(
  Famine_start ~ GVF:i(Region) + nino34 | Region + Decade,
  data     = famines,
  panel.id = c("Region", "Year"),
  vcov     = DK ~ Year
)

# Extract GVF × Region coefficients
b_vr   <- coef(m_volc_region)
V_vr   <- vcov(m_volc_region)
idx_vr <- grep("GVF", names(b_vr))

coef_vr <- tibble(
  term     = names(b_vr)[idx_vr],
  estimate = b_vr[idx_vr],
  se       = sqrt(diag(V_vr)[idx_vr])
) %>%
  mutate(
    lo95 = estimate - 1.960 * se,
    hi95 = estimate + 1.960 * se,
    lo90 = estimate - 1.645 * se,
    hi90 = estimate + 1.645 * se,
    Region = term %>%
      gsub("GVF:i\\(Region\\)::", "", .) %>%
      gsub("GVF:Region::",        "", .) %>%
      trimws()
  )

print(coef_vr[, c("Region", "estimate", "lo95", "hi95")])

panel_e <- ggplot(coef_vr, aes(x = Region, y = estimate)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35,
             linetype = "dashed") +
  geom_errorbar(aes(ymin = lo95, ymax = hi95),
                width = 0.0, linewidth = 0.55, colour = "tomato") +
  geom_errorbar(aes(ymin = lo90, ymax = hi90),
                width = 0.0, linewidth = 1.5, colour = "tomato") +
  geom_point(size = 2.4, colour = "tomato", shape = 17) +
  annotate("rect",
           xmin = 0.5, xmax = 1.5,
           ymin = -Inf, ymax = Inf,
           fill = "tomato", alpha = 0.07) +
  labs(x = NULL, y = "V.F. effect on famine onset prob.",
       title = "V.F. effect by region (controlling for ENSO)") +
  nature_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10))

panel_e

fig2 <- (panel_d | panel_e) +
  plot_layout(widths = c(1.4, 1.2)) +
  plot_annotation(tag_levels = list(c("d", "e"))) &
  theme(plot.tag   = element_text(face = "bold", size = 20, colour = "black"),
        plot.title = element_text(face = "bold", size = 13))

fig2

out_dir <- file.path(ROOT_DIR, "analysis/output/figures/extended data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(out_dir, "figED_volcanic.pdf"),
  plot   = fig,
  width  = 12,
  height = 9,
  units  = "in",
  device = cairo_pdf
)

message("Saved → figED_volcanic.pdf")
