###############################################################################
# Fed-GEE VISN application - forest plot, Stat Med version (2026-09-29)
#
# Three analyses of the same K = 18 VISN data:
#   Pooled GEE        - all data in one place, VISN-clustered sandwich, z.
#                       (Identical to uncorrected Fed-GEE: same beta, same SE.)
#   Fed-GEE (ours)    - Kauermann-Carroll score-space correction + t with
#                       score-space Bell-McCaffrey df. The recommended analysis.
#   Meta-analytic GEE - each VISN fits its own GEE; fixed-effect inverse-
#                       variance pool, z. Uses only the VISNs that pass the
#                       per-VISN fitting screen.
# Pooled and Fed-GEE share point estimates; meta-GEE does not.
# Source: results/ (VINCI run 2026-09-29 of visn_analysis.R).
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

# Run from inst/application:   Rscript visn_forest.R
IN_DIR <- "results"
OUT <- "../manuscript/figures/visn_forest"

x <- readRDS(file.path(IN_DIR, "2026_09_27_visn_results_central.RDS"))
est <- x$centralized
meta_fit <- read.csv(file.path(IN_DIR, "2026_09_27_visn_meta.csv"))

lab_pooled <- "Pooled GEE"
lab_fed <- "Fed-GEE: KC + Bell-McCaffrey df (recommended)"
lab_meta <- "Meta-analytic GEE"
keep <- c(lab_fed, lab_pooled, lab_meta)

d <- bind_rows(
  est %>% filter(variant == "site", inference == "z") %>%
    transmute(term, OR, CI_low, CI_high, Model = lab_pooled),
  est %>% filter(variant == "kc", inference == "t_bm") %>%
    transmute(term, OR, CI_low, CI_high, Model = lab_fed),
  meta_fit %>% transmute(term, OR, CI_low, CI_high, Model = lab_meta)
) %>% mutate(Model = factor(Model, levels = rev(keep))) # last level plots on top

## validated categorical palette (dataviz defaults, light mode: all checks
## pass; aqua < 3:1 contrast -> legend + shape carry identity too).
pal <- setNames(c("#eb6834", "#2a78d6", "#1baf7a"), keep)
shp <- setNames(c(16, 16, 17), keep)

## ---- term -> variable / level / category ---------------------------------
meta <- tribble(
  ~term, ~Variable, ~Level, ~Category,
  "female1", "Female\n(vs male)", NA, "Demographics",
  "age651", "Age 65+\n(vs 18-64)", NA, "Demographics",
  "auditc_low1", "AUDIT-C no or low use\n(vs unhealthy or higher)", NA, "Alcohol severity",
  "oud1", "OUD diagnosis\n(vs none)", NA, "Alcohol severity",
  "can_high1", "CAN score >= 60th pctile\n(top two quintiles vs lower)", NA, "Medical complexity",
  "frail1", "Frailty, mild or worse\n(vs nonfrail/prefrail)", NA, "Medical complexity",
  "icu_any1", "ICU stay\n(vs none)", NA, "Medical complexity",
  "psych1", "Psychiatry service\n(vs med/surg)", NA, "Hospital process",
  "consult31", "Addiction consult\n(ref: none)", "Clinician", "Hospital process",
  "consult32", "Addiction consult\n(ref: none)", "Non-clinician", "Hospital process",
  "ama1", "Discharge against advice\n(vs routine)", NA, "Hospital process"
)
d <- d %>% inner_join(meta, by = "term")
d$row_label <- ifelse(is.na(d$Level), d$Variable, paste0("    ", d$Level))

## header rows for the one multi-level variable (addiction consult)
headers <- meta %>%
  filter(!is.na(Level)) %>%
  distinct(Variable, Category) %>%
  transmute(
    term = paste0("__h_", Variable), Variable, Level = NA_character_,
    Category, row_label = Variable,
    Model = factor(NA, levels = rev(keep)),
    OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_
  )
plot_df <- bind_rows(d, headers)

## ---- ordering -------------------------------------------------------------
cat_levels <- c("Demographics", "Alcohol severity", "Medical complexity", "Hospital process")
plot_df$Category <- factor(plot_df$Category, levels = cat_levels)
row_seq <- c()
seen <- character(0)
for (i in seq_len(nrow(meta))) {
  v <- meta$Variable[i]
  lv <- meta$Level[i]
  if (!is.na(lv) && !(v %in% seen)) {
    row_seq <- c(row_seq, v)
    seen <- c(seen, v)
  }
  row_seq <- c(row_seq, if (is.na(lv)) v else paste0("    ", lv))
}
plot_df$row_label <- factor(plot_df$row_label, levels = rev(unique(row_seq)))
plot_df$is_header <- is.na(plot_df$Model)

## direct labels: OR on the recommended (KC) variant only, to avoid clutter
lab_df <- plot_df %>%
  filter(Model == lab_fed) %>%
  mutate(txt = sprintf("%.2f", OR))

## ---- caption --------------------------------------------------------------
N <- x$diagnostics$n_patients # unique veterans (GEE cluster count)
Kv <- x$diagnostics$n_sites # VISNs
# Hospitalization (row) count is NOT stored in the RDS. Set it from the VINCI
# run: N_HOSP <- nrow(df) inside main(). AoIM Table 1 reports 29,041; confirm
# this run matches before the figure is shared.
N_HOSP <- format(29041, big.mark = ",")
cap_raw <- paste(
  sprintf(
    "Note.\n1. Marginal logistic GEE of in-hospital MAUD initiation across %s hospitalizations of %s unique veterans in %d VISNs (the unit of federation). Pooled GEE and Fed-GEE give identical odds ratios; Fed-GEE differs only in its interval.",
    N_HOSP, format(N, big.mark = ","), Kv
  ),
  "2. Pooled GEE: VISN-clustered sandwich with a normal reference distribution. Fed-GEE: Kauermann-Carroll score-space correction with a t-distribution reference on Bell-McCaffrey degrees of freedom, computed from the VISN summaries alone. Meta-analytic GEE: each VISN fits its own GEE and the estimates are inverse-variance pooled.",
  sep = "\n"
)
# Notes are drawn below the plot as a justified text block (ggplot captions
# cannot justify). Numbered notes get a hanging indent.
notes <- strsplit(cap_raw, "\n")[[1]]

# Lay out justified lines; measures text on the device that is drawing, so
# the PDF and PNG each get exact spacing.
layout_notes <- function(paras, width_in, gp, lineheight = 1.3) {
  sw <- function(t) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(t, gp = gp)),
      "in",
      valueOnly = TRUE
    )
  }
  space <- sw("a a") - sw("aa")
  line_h <- gp$fontsize * lineheight / 72
  grobs <- list()
  y <- 0
  for (para in paras) {
    num <- regmatches(para, regexpr("^[0-9]+\\. ", para))
    body <- sub("^[0-9]+\\. ", "", para)
    indent <- if (length(num)) sw("0. ") + 0.02 else 0
    if (length(num)) {
      grobs[[length(grobs) + 1]] <- grid::textGrob(
        trimws(num),
        x = grid::unit(0, "in"), y = grid::unit(1, "npc") - grid::unit(y, "in"),
        hjust = 0, vjust = 1, gp = gp
      )
    }
    words <- strsplit(body, " +")[[1]]
    wl <- vapply(words, sw, numeric(1))
    avail <- width_in - indent
    lines <- list()
    cur <- integer(0)
    cur_w <- 0
    for (i in seq_along(words)) {
      add <- if (length(cur)) space + wl[i] else wl[i]
      if (length(cur) && cur_w + add > avail) {
        lines[[length(lines) + 1]] <- cur
        cur <- i
        cur_w <- wl[i]
      } else {
        cur <- c(cur, i)
        cur_w <- cur_w + add
      }
    }
    lines[[length(lines) + 1]] <- cur
    for (li in seq_along(lines)) {
      idx <- lines[[li]]
      last <- li == length(lines) || length(idx) == 1
      gap <- if (last) space else (avail - sum(wl[idx])) / (length(idx) - 1)
      xs <- indent + c(0, cumsum(wl[idx] + gap))[seq_along(idx)]
      grobs[[length(grobs) + 1]] <- grid::textGrob(
        words[idx],
        x = grid::unit(xs, "in"),
        y = grid::unit(1, "npc") - grid::unit(y, "in"),
        hjust = 0, vjust = 1, gp = gp
      )
      y <- y + line_h
    }
    y <- y + line_h * 0.25
  }
  list(grobs = grobs, height_in = y)
}

justified_notes_grob <- function(paras, width_in, fontsize = 8,
                                 family = "Helvetica", col = "#1a1a1a") {
  gp <- grid::gpar(fontsize = fontsize, fontfamily = family, col = col)
  h <- layout_notes(paras, width_in, gp)$height_in
  g <- grid::gTree(
    paras = paras, width_in = width_in, gp_notes = gp,
    cl = "justnotes"
  )
  list(grob = g, height_in = h)
}
makeContent.justnotes <- function(x) {
  grid::setChildren(x, do.call(
    grid::gList,
    layout_notes(x$paras, x$width_in, x$gp_notes)$grobs
  ))
}

## ---- plot -----------------------------------------------------------------
p <- ggplot(plot_df, aes(x = OR, y = row_label, color = Model)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = .4, color = "#b0b0b0") +
  geom_linerange(aes(xmin = CI_low, xmax = CI_high),
    position = position_dodge(width = .62),
    linewidth = .7, na.rm = TRUE
  ) +
  geom_point(aes(shape = Model),
    position = position_dodge(width = .62),
    size = 2.1, na.rm = TRUE
  ) +
  geom_text(
    data = lab_df, aes(x = CI_high, label = txt),
    hjust = -0.25, size = 2.7, color = "#3a3a3a", na.rm = TRUE
  ) +
  scale_color_manual(values = pal, breaks = names(pal)) +
  scale_shape_manual(values = shp, breaks = names(shp)) +
  scale_x_log10(
    breaks = c(.6, .75, 1, 1.5, 2, 3),
    labels = c("0.6", "0.75", "1", "1.5", "2", "3"),
    expand = expansion(mult = c(.02, .10))
  ) +
  facet_grid(Category ~ .,
    scales = "free_y", space = "free_y",
    labeller = labeller(Category = label_wrap_gen(10))
  ) +
  labs(
    x = "Adjusted odds ratio (95% CI) of MAUD initiation",
    y = NULL, color = NULL, shape = NULL,
    title = "MAUD initiation at K = 18 VISNs: pooled, federated, and meta-analytic GEE",
  ) +
  theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid.major.x = element_line(color = "#ececec", linewidth = .4),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(10, "pt"),
    strip.text.y = element_text(
      angle = 0, hjust = 0.5, face = "bold",
      size = 8.5, color = "white"
    ),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text.y = element_text(color = "#222222", size = 8.5, hjust = 0),
    axis.text.x = element_text(color = "#555555", size = 8.5),
    axis.title.x = element_text(size = 9, color = "#333333", margin = margin(t = 6)),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 8),
    legend.key.spacing.x = unit(6, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(
      size = 11.5, face = "bold", color = "#1a1a1a",
      margin = margin(b = 2)
    ),
    plot.title.position = "plot",
    plot.caption = element_text(
      hjust = 0, size = 8, color = "#1a1a1a",
      lineheight = 1.25, margin = margin(t = 12)
    ),
    plot.caption.position = "plot",
    plot.margin = margin(12, 16, 10, 12)
  ) +
  guides(shape = "none", color = guide_legend(
    nrow = 1, reverse = FALSE,
    override.aes = list(linewidth = 1.4, size = 2.6, shape = unname(shp))
  ))

FIG_W <- 8.5
FIG_H <- 11
MARG_L <- 12 / 72
MARG_R <- 16 / 72
nt <- justified_notes_grob(notes, width_in = FIG_W - MARG_L - MARG_R)
notes_panel <- grid::gTree(
  children = grid::gList(nt$grob),
  vp = grid::viewport(
    x = grid::unit(MARG_L, "in"), y = grid::unit(1, "npc") - grid::unit(4 / 72, "in"),
    just = c(0, 1), width = grid::unit(FIG_W - MARG_L - MARG_R, "in"),
    height = grid::unit(nt$height_in, "in")
  )
)
full <- gridExtra::arrangeGrob(ggplotGrob(p), notes_panel,
  ncol = 1,
  heights = grid::unit.c(grid::unit(1, "null"), grid::unit(nt$height_in + 0.25, "in"))
)

ggsave(paste0(OUT, "_standalone.pdf"), full, width = FIG_W, height = FIG_H, bg = "white")
ggsave(paste0(OUT, "_standalone.png"), full, width = FIG_W, height = FIG_H, dpi = 200, bg = "white")
cat("saved", paste0(OUT, "_standalone.pdf/.png"), "\n")

## ---- manuscript version: no title, no notes (they go in the LaTeX caption) ----
p_paper <- p + labs(title = NULL) +
  theme(plot.margin = margin(4, 8, 4, 4))
ggsave(paste0(OUT, "_paper.pdf"), p_paper, width = 7, height = 7.6, bg = "white")
ggsave(paste0(OUT, "_paper.png"), p_paper, width = 7, height = 7.6, dpi = 200, bg = "white")
cat("saved", paste0(OUT, "_paper.pdf/.png"), "\n")
