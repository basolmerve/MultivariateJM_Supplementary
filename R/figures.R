library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

pd_data_long <- read_delim(file = "data/longitudinal_data84.txt", delim = "\t", col_names = TRUE, 
                           locale = locale(decimal_mark = "."))

pd_data_additional <- read_delim("data/additional_data.csv", delim = ",", locale = locale(decimal_mark = "."))

pd_data <- pd_data_long %>%
  full_join(pd_data_additional) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(preHD = factor(preHD, levels = c(0, 1), labels = c("No", "Yes")),
         highlow = factor(highlow, levels = c(1, 2), labels = c("Low", "High"))) %>%
  as.data.frame

## longitudinal ##
pd_data_train <- pd_data %>%
  filter(complete.cases(time, event, gender, age, BKI, preHD, ill_count, peritonitrate, ALB, 
                        BUN, KR, CA, P)) %>% 
  filter(!(id %in% c(13, 4)))  # 2 gozlem veri setinden cikarildi.

# Counting process
counting <- function(times, survtime){
  c(times[-1], survtime)
}

pd_data_train <- pd_data_train %>%
  group_by(id) %>%
  mutate(tstart = time, tstop = counting(time, unique(surv)),
         statusC = if_else(
           event == 0, 
           rep(0, length(time)),
           rep(c(0, 1), c(length(time) - 1, 1))
         ))

test1 <- pd_data_long[pd_data_long$id == 4, ]
test2 <- pd_data_long[pd_data_long$id == 13, ]
TestSamples <- rbind(test1, test2)

#### Graphs #####
point_size <- 1
line_width <- .5
jitter_amount <- 4
span_amount <- .75
default_theme <- theme_bw(base_size = 8) + 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 2.5, b = 2.5)),
    axis.text.y = element_text(margin = margin(r = 2.5, l = 2.5))
  )

alb <- pd_data_train %>%
  dplyr::select(one_of(c("id", "event2", "time", "ALB"))) %>%
  ggplot(aes(x = time, y = ALB)) + 
  default_theme +
  geom_jitter(width = jitter_amount, size = point_size, colour = rgb(0, 0, 0, 0.1)) + 
  geom_smooth(method = "loess", se = FALSE, colour = "red", lwd = line_width, span = span_amount) +
  geom_smooth(data = TestSamples, aes(x = time, y = ALB), method = "loess", se = FALSE, colour = "green2", lwd = line_width, span = span_amount) +
  facet_wrap(~ event2) + 
  labs(x = "Time (month)", y = "Albumin (mg/dL)")

bun <- pd_data_train %>%
  dplyr::select(one_of(c("id", "event2", "time", "BUN"))) %>%
  ggplot(aes(x = time, y = BUN)) + 
  default_theme +
  geom_jitter(width = jitter_amount, size = point_size, colour = rgb(0, 0, 0, 0.1)) + 
  geom_smooth(method = "loess", se = FALSE, colour = "red", lwd = line_width, span = span_amount) +
  geom_smooth(data = TestSamples, aes(x = time, y = BUN), method = "loess", se = FALSE, colour = "green2", lwd = line_width, span = span_amount) +
  facet_wrap(~ event2) + 
  labs(x = "Time (month)", y = "BUN (mg/dL)")

kr <- pd_data_train %>%
  dplyr::select(one_of(c("id", "event2", "time", "KR"))) %>%
  ggplot(aes(x = time, y = KR)) + 
  default_theme + 
  geom_jitter(width = jitter_amount, size = point_size, colour = rgb(0, 0, 0, 0.1)) + 
  geom_smooth(method = "loess", se = FALSE, colour = "red", lwd = line_width, span = span_amount) +
  geom_smooth(data = TestSamples, aes(x = time, y = KR), method = "loess", se = FALSE, colour = "green2", lwd = line_width, span = span_amount) +
  facet_wrap(~ event2) + 
  labs(x = "Time (month)", y = "Creatinine (mg/dL)")


ca <- pd_data_train %>%
  dplyr::select(one_of(c("id", "event2", "time", "CA"))) %>%
  ggplot(aes(x = time, y = CA)) + 
  default_theme + 
  geom_jitter(width = jitter_amount, size = point_size, colour = rgb(0, 0, 0, 0.1)) + 
  geom_smooth(method = "loess", se = FALSE, colour = "red", lwd = line_width, span = span_amount) +
  geom_smooth(data = TestSamples, aes(x = time, y = CA), method = "loess", se = FALSE, colour = "green2", lwd = line_width, span = span_amount) +
  facet_wrap(~ event2) + 
  labs(x = "Time (month)", y = "Calcium (mg/dL)")

p <- pd_data_train %>%
  dplyr::select(one_of(c("id", "event2", "time", "P"))) %>%
  ggplot(aes(x = time, y = P)) + 
  default_theme + 
  geom_jitter(width = jitter_amount, size = point_size, colour = rgb(0, 0, 0, 0.1)) + 
  geom_smooth(method = "loess", se = FALSE, colour = "red", lwd = line_width, span = span_amount) +
  geom_smooth(data = TestSamples, aes(x = time, y = P), method = "loess", se = FALSE, colour = "green2", lwd = line_width, span = span_amount) +
  facet_wrap(~ event2) + 
  labs(x = "Time (month)", y = "Phosphate (mg/dL)")

figure <- ggarrange(alb, bun, kr, ca, p,
                    ncol = 3, nrow = 2)

ggsave(filename = "figure/fig1.eps", plot = figure, device = cairo_ps, width = 7.5, height = 4)

## longitudinal ##
longitudinalData <- pd_data_long %>% 
  filter(id %in% c(4, 13)) %>%
  dplyr::select(one_of(c("id", "time", "ALB", "BUN", "KR", "CA", "P"))) %>%
  mutate(id = as.factor(id)) %>%
  as_tibble() 

map <- function(x, min.y = 0, max.y = 1, xmin = NULL, xmax = NULL, ..){
  if (is.null(xmin)){
    xmin <- min(x)
  }
  
  if (is.null(xmax)){
    xmax <- max(x)
  }
  
  ratio <- (x - xmin) / (xmax - xmin)
  min.y + ratio * (max.y - min.y)
}

# prob4 <- read.table("ID4.txt", header = T, sep = "\t")
# prob13 <- read.table("ID13.txt", header = T, sep = "\t")
prob4 <- read.csv(file = "data/ID4.csv")
prob13 <- read.csv(file = "data/ID13.csv")

survProbData <- bind_rows(prob4, prob13) %>%
  mutate(id = as.factor(id),
         method = as.factor(method)) %>%
  as_tibble()

#### IDs 4 and 13 
thresholdALB <- bind_rows(
  tibble(id = 4, x = c(0, 48), y = c(3.1, 3.1)),
  tibble(id = 13, x = c(0, 42), y = c(3.62, 3.62))
) %>%
  mutate(id = as.factor(id))

thresholdBUN <- bind_rows(
  tibble(id = 4, x = c(0, 48), y = c(49.22, 49.22)),
  tibble(id = 13, x = c(0, 42), y = c(64.41, 64.41))
) %>%
  mutate(id = as.factor(id))

thresholdKR <- bind_rows(
  tibble(id = 4, x = c(0, 48), y = c(7.14, 7.14)),
  tibble(id = 13, x = c(0,42), y = c(6.66, 6.66))
) %>%
  mutate(id = as.factor(id))

thresholdCA <- bind_rows(
  tibble(id = 4, x = c(0, 48), y = c(8.91, 8.91)),
  tibble(id = 13, x = c(0, 42), y = c(8.61, 8.61))
) %>%
  mutate(id = as.factor(id))

thresholdP <- bind_rows(
  tibble(id = 4, x = c(0, 48), y = c(3.44, 3.44)),
  tibble(id = 13, x = c(0, 42), y = c(5.5, 5.5))
) %>%
  mutate(id = as.factor(id))


######  FIGURE 2 ########
vlineData <- longitudinalData %>%
  group_by(id) %>%
  summarise(vline = max(time))

removeStrip <- theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
)

facetLabels <- c("Patient A", "Patient B")
names(facetLabels) <- c("4", "13")

# plot options
baseFontSize <- 8
pointSize <- 1
lineWidth <- .5
sec.axis.title <- "Survival probs"

default_theme <- theme_bw(base_size = baseFontSize)  +
  theme(panel.grid = element_blank(), legend.position = "top", 
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(margin = margin(l = 2.5, r = 2.5)),
        axis.text.y.right = element_text(margin = margin(l = 2.5, r = 2.5)),
        axis.title.y.right = element_text(margin = margin(l = 5)),
        axis.title.x = element_text(margin = margin(t = 5)))

plotALB <- ggplot(longitudinalData, aes(x = time, y = ALB)) + 
  default_theme +
  geom_point(size = pointSize, colour = rgb(0, 0, 0, 0.6)) +
  geom_line(linewidth = lineWidth) + 
  geom_line(data = thresholdALB, 
            mapping = aes(x = x, y = y), colour = "red", linetype = 2, linewidth = lineWidth) +
  geom_line(data = survProbData, 
            mapping = aes(x = t, y = map(ALB, min.y = 0, max.y = 5.5, xmin = 0, xmax = 1), 
                          colour = method), 
            linewidth = lineWidth) + 
  geom_vline(aes(xintercept = vline), data = vlineData, lty = 2, colour = "gray30") +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 130, 20)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks = map(seq(0, 1, 0.2), min.y = 0, 
                                                      max.y = 5.5, xmin = 0, xmax = 1),
                                         labels = formatC(seq(0, 1, 0.2), digits = 1, format = "f"), 
                                         name = sec.axis.title), limits = c(0, 5.5), breaks = 0:5) + 
  labs(y = "Albumin (mg/dL)", x = NULL) + 
  facet_grid(cols = vars(id), labeller = labeller(id = facetLabels)) + 
  guides(colour = guide_legend(title = "Method")) + 
  removeStrip


plotBUN <- ggplot(longitudinalData, aes(x = time, y = BUN)) + 
  default_theme + 
  geom_point(size = pointSize, colour = rgb(0, 0, 0, 0.6)) +
  geom_line(linewidth = lineWidth) + 
  geom_line(data = thresholdBUN, 
            mapping = aes(x = x, y = y), colour = "red", linetype = 2, linewidth = lineWidth) +
  geom_line(data = survProbData, 
            mapping = aes(x = t, y = map(BUN, min.y = 20, max.y = 95, xmin = 0, xmax = 1), colour = method), 
            linewidth = lineWidth) + 
  geom_vline(aes(xintercept = vline), data = vlineData, lty = 2, colour = "gray30") +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 130, 20)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks = map(seq(0, 1, 0.2), min.y = 20, max.y = 95, xmin = 0, xmax = 1),
                                         labels = formatC(seq(0, 1, 0.2), digits = 1, format = "f"), 
                                         name = sec.axis.title), limits = c(20, 95)) + 
  labs(y = "BUN (mg/dL)", x = NULL) + 
  facet_grid(cols = vars(id)) + removeStrip +
  guides(colour = "none")


plotKR <- ggplot(longitudinalData, aes(x = time, y = KR)) + 
  default_theme + 
  geom_point(size = pointSize, colour = rgb(0, 0, 0, 0.6)) + 
  geom_line(linewidth = lineWidth) + 
  geom_line(data = thresholdKR, 
            mapping = aes(x = x, y = y), colour = "red", linetype = 2, linewidth = lineWidth) +
  geom_line(data = survProbData, 
            mapping = aes(x = t, y = map(KR, min.y = 0, max.y = 18, xmin = 0, xmax = 1), colour = method), 
            linewidth = lineWidth) + 
  geom_vline(aes(xintercept = vline), data = vlineData, lty = 2, colour = "gray30") +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 130, 20)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks = map(seq(0, 1, 0.2), min.y = 0, max.y = 18, xmin = 0, xmax = 1),
                                         labels = formatC(seq(0, 1, 0.2), digits = 1, format = "f"), 
                                         name = sec.axis.title), limits = c(0, 18)) + 
  labs(y = "Creatinine (mg/dL)", x = NULL) + 
  facet_grid(cols = vars(id)) + removeStrip +
  guides(colour = "none")


plotCA <- ggplot(longitudinalData, aes(x = time, y = CA)) + 
  default_theme + 
  geom_point(size = pointSize, colour = rgb(0, 0, 0, 0.6)) + 
  geom_line(linewidth = lineWidth) + 
  geom_line(data = thresholdCA, 
            mapping = aes(x = x, y = y), colour = "red", linetype = 2, linewidth = lineWidth) +
  geom_line(data = survProbData, 
            mapping = aes(x = t, y = map(CA, min.y = 3, max.y = 12, xmin = 0, xmax = 1), colour = method), 
            linewidth = lineWidth) + 
  geom_vline(aes(xintercept = vline), data = vlineData, lty = 2, colour = "gray30") +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 130, 20)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks = map(seq(0, 1, 0.2), min.y = 3, max.y = 12, xmin = 0, xmax = 1),
                                         labels = formatC(seq(0, 1, 0.2), digits = 1, format = "f"), 
                                         name = sec.axis.title), limits = c(3, 12)) + 
  labs(y = "Calcium (mg/dL)", x = NULL) + 
  facet_grid(cols = vars(id)) + removeStrip +
  guides(colour = "none")


plotP <- ggplot(longitudinalData, aes(x = time, y = P)) + 
  default_theme + 
  geom_point(size = pointSize, colour = rgb(0, 0, 0, 0.6)) + 
  geom_line(linewidth = lineWidth) + 
  geom_line(data = thresholdP, 
            mapping = aes(x = x, y = y), colour = "red", linetype = 2, linewidth = lineWidth) +
  geom_line(data = survProbData, 
            mapping = aes(x = t, y = map(P, min.y = 0, max.y = 12, xmin = 0, xmax = 1), colour = method), 
            linewidth = lineWidth) + 
  geom_vline(aes(xintercept = vline), data = vlineData, lty = 2, colour = "gray30") +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 130, 20)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks = map(seq(0, 1, 0.2), min.y = 0, max.y = 12, xmin = 0, xmax = 1),
                                         labels = formatC(seq(0, 1, 0.2), digits = 1, format = "f"), 
                                         name = sec.axis.title), limits = c(0, 12)) + 
  labs(y = "Phosphate (mg/dL)", x = "Time (months)") + 
  facet_grid(cols = vars(id)) + removeStrip +
  guides(colour = "none")


allPlots <- ggarrange(plotALB, plotBUN, plotKR, plotCA, plotP,
                      ncol = 1, nrow = 5)

ggsave(plot = allPlots, filename = "figure/dynPred.eps", device = cairo_ps, height = 8, width = 5)

############ FIGURE 3 - AUC GRAPHS ###################
AUC_data <- read.table(file = "data/AUC_data.txt",header = T,sep = "\t")

point_size <- 1
line_width <- .5
jitter_amount <- 4
span_amount <- .75
default_theme <- theme_bw(base_size = 8) + 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 2.5, b = 2.5)),
    axis.text.y = element_text(margin = margin(r = 2.5, l = 2.5))
  )

AUC_ALB_graph <- ggplot(AUC_data, aes(x = time, y = ALB, group = Method, colour = Method)) +
  geom_line(linewidth = line_width) + 
  geom_line(aes(y = 0.5), linetype = "longdash", linewidth = line_width, colour = "black") +
  geom_point(size = point_size) + 
  default_theme +
  scale_x_continuous(breaks = seq(0, 66, 6)) + 
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.5, 1, 0.1)) + 
  labs(y = "Area Under Curve (AUC)", x = "Time (months)") + 
  guides(colour = guide_legend(title = "Method")) + 
  ggtitle("Serum Albumin (mg/dL)")


AUC_BUN_graph <- ggplot(AUC_data, aes(x = time, y = BUN, group = Method, colour = Method)) +
  geom_line(linewidth = line_width) + 
  geom_line(aes(y = 0.5), linetype = "longdash", linewidth = line_width, colour = "black") +
  geom_point(size = point_size) + 
  default_theme + 
  scale_x_continuous(breaks = seq(0, 66, 6)) + 
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.5, 1, 0.1)) + 
  labs(y = "Area Under Curve (AUC)", x = "Time (months)") + 
  guides(colour = guide_legend(title = "Method")) + 
  ggtitle("BUN (mg/dL)")

AUC_KR_graph <- ggplot(AUC_data, aes(x = time, y = KR, group = Method, colour = Method)) +
  geom_line(linewidth = line_width) + 
  geom_line(aes(y = 0.5), linetype = "longdash", linewidth = line_width, colour = "black") +
  geom_point(size = point_size) + 
  default_theme + 
  scale_x_continuous(breaks = seq(0, 66, 6)) + 
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.5, 1, 0.1)) + 
  labs(y = "Area Under Curve (AUC)", x = "Time (months)") + 
  guides(colour = guide_legend(title = "Method")) + 
  ggtitle("Creatinine (mg/dL)")


AUC_CA_graph <- ggplot(AUC_data, aes(x = time, y = CA, group = Method, colour = Method)) +
  geom_line(linewidth = line_width) + 
  geom_line(aes(y = 0.5), linetype = "longdash", linewidth = line_width, colour = "black") +
  geom_point(size = point_size) + 
  default_theme + 
  scale_x_continuous(breaks = seq(0, 66, 6)) + 
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.5, 1, 0.1)) + 
  labs(y = "Area Under Curve (AUC)", x = "Time (months)") + 
  guides(colour = guide_legend(title = "Method")) + 
  ggtitle("Calcium (mg/dL)")

AUC_P_graph <- ggplot(AUC_data, aes(x = time, y = P, group = Method, colour = Method)) +
  geom_line(linewidth = line_width) + 
  geom_line(aes(y = 0.5), linetype = "longdash", linewidth = line_width, colour = "black") +
  geom_point(size = point_size) + 
  default_theme + 
  scale_x_continuous(breaks = seq(0, 66, 6)) + 
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.5, 1, 0.1)) + 
  labs(y = "Area Under Curve (AUC)", x = "Time (months)") + 
  guides(colour = guide_legend(title = "Method")) + 
  ggtitle("Phosphate (mg/dL)")


AUCPlots <- ggarrange(AUC_ALB_graph, AUC_BUN_graph, AUC_KR_graph, AUC_CA_graph,
                      AUC_P_graph, ncol = 2, nrow = 3, legend = "top", common.legend = TRUE)

ggsave(plot = AUCPlots, filename = "figure/AUCPlots.eps", device = cairo_ps, height = 8, width = 5)
