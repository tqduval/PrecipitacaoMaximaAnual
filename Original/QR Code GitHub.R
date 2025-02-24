
# QRCODE ------------------------------------------------------------------

install.packages("qrcode")
library(qrcode)

generate_svg(qrcode = qr_code("https://github.com/tqduval/PDMaxBrasil"),
             filename = "Dados Gerados/git.qrcode.svg")

plot(qr_code("https://github.com/tqduval/PDMaxBrasil"))