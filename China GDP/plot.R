library(readxl)
library(ggplot2)

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

D <- read_excel("service_ratio_plot.xlsx")
D$Year <- seq.Date(as.Date("2000/1/1"), as.Date("2007/1/1"), "year")

name_vec <- colnames(D)[-1]
col_vec <- gg_color_hue(5)
df_plot <- reshape2::melt(D, id = "Year")

pic <- ggplot(data = df_plot) + 
    geom_line(mapping = aes(
        x = Year,
        y = value * 100,
        color = variable,
        size = variable
    )) +
    scale_colour_manual(breaks = name_vec,
                        values = col_vec,
                        labels = name_vec) +
    scale_size_manual(breaks = name_vec,
                      values = c(1.1, rep(0.75, 4)),
                      labels = name_vec) +
    scale_x_date(breaks = D$Year, labels = as.character(lubridate::year(D$Year))) +
    labs(x = "Year", y = "Tertiary Industry (%)") +
    theme(
        panel.background =  element_blank(),
        panel.border = element_rect(
            linetype = 1,
            colour = "black",
            fill = NA
        ),
        panel.grid.major = element_line(linetype = 2, color = "grey90"),
        legend.title = element_blank(),
        legend.position = "bottom"
    )


