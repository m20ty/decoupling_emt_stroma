# g <- ggplot(
#     data.frame(
#         pos = 1:length(labels_vec),
#         label = labels_vec
#     )
# ) + theme(
#     panel.background = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks.length = unit(0, 'pt'),
#     axis.ticks = element_blank(),
#     axis.title = element_text(size = axis_title_size)
# )

# ggpubr::geom_bracket(
#     xmin = 0,
#     xmax = 20,
#     y.position = 0,
#     label = 'pEMT',
#     vjust = 2
# ) +
# ggpubr::geom_bracket(
#     xmin = 80,
#     xmax = 100,
#     y.position = 0,
#     label = 'CAF',
#     vjust = 2
# )

ggplot(
    data.frame(x = 1:100, y = seq(-1, 1, length.out = 100))
) +
    theme(
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt')
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) +
    geom_segment(aes(x = 0, xend = 20, y = 0, yend = 0)) +
    geom_segment(aes(x = 80, xend = 100, y = 0, yend = 0)) +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1)) +
    geom_segment(aes(x = 20, xend = 20, y = 0, yend = 1)) +
    geom_segment(aes(x = 80, xend = 80, y = 0, yend = 1)) +
    geom_segment(aes(x = 100, xend = 100, y = 0, yend = 1)) +
    annotate(geom = 'text', x = 10, y = 0, label = 'pEMT', vjust = 1.5) +
    annotate(geom = 'text', x = 90, y = 0, label = 'CAF', vjust = 1.5)

pemt_caf_brackets <- function(
    left_bracket_xmax = 20,
    right_bracket_xmin = 80,
    brackets_y = 0,
    tips_yend = 1,
    labels = c('pEMT', 'CAF'),
    labels_vjust = 1.5
) {
    
    ggplot(
        data.frame(x = 1:100, y = seq(-1, 1, length.out = 100))
    ) +
        theme(
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, 'pt')
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) +
        geom_segment(
            aes(
                x = 0,
                xend = left_bracket_xmax,
                y = brackets_y,
                yend = brackets_y
            )
        ) +
        geom_segment(
            aes(
                x = right_bracket_xmin,
                xend = 100,
                y = brackets_y,
                yend = brackets_y
            )
        ) +
        geom_segment(
            aes(
                x = 0,
                xend = 0,
                y = brackets_y,
                yend = tips_yend
            )
        ) +
        geom_segment(
            aes(
                x = left_bracket_xmax,
                xend = left_bracket_xmax,
                y = brackets_y,
                yend = tips_yend
            )
        ) +
        geom_segment(
            aes(
                x = right_bracket_xmin,
                xend = right_bracket_xmin,
                y = brackets_y,
                yend = tips_yend
            )
        ) +
        geom_segment(
            aes(
                x = 100,
                xend = 100,
                y = brackets_y,
                yend = tips_yend
            )
        ) +
        annotate(
            geom = 'text',
            x = left_bracket_xmax/2,
            y = brackets_y,
            label = labels[1],
            vjust = labels_vjust
        ) +
        annotate(
            geom = 'text',
            x = right_bracket_xmin + (100 - right_bracket_xmin)/2,
            y = brackets_y,
            label = labels[2],
            vjust = labels_vjust
        )
    
}
