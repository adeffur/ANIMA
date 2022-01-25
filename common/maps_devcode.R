devtools::install_github("jmsallan/BAdatasetsSpatial")
library(BAdatasetsSpatial)
library(tidyverse)
names(WorldMap1_110)
ggplot(WorldMap1_110) + 
  geom_sf() +
  theme_void()
WorldMap1_10 %>%
  filter(ISO_A3 == "ESP") %>%
  ggplot + 
  geom_sf() +
  theme_void()
WorldMap1_10 %>%
  filter(ISO_A3 == "GER") %>%
  ggplot + 
  geom_sf() +
  theme_void()
WorldMap1_110 %>%
  mutate(map_color = as.factor(MAPCOLOR7)) %>%
  ggplot + 
  geom_sf(aes(fill= map_color)) +
  theme_void() +
  scale_fill_brewer(palette = "Set1") +
  theme(panel.background = element_rect(fill = "#CCE5FF"),
        legend.position = "none")
gdp <- WB_gdp_per_capita %>%
  filter(year == 2019) %>%
  mutate(gdp = cut(value, breaks = 12))
gdp %>%
  group_by(gdp) %>%
  summarise(n = n(), .groups = "drop")
gdp <- gdp %>%
  mutate(gdp = fct_lump(gdp, 5))
gdp %>%
  group_by(gdp) %>%
  summarise(n = n(), .groups = "drop")
WorldMap1_110 %>%
  ggplot +
  geom_sf(aes(fill = gdp)) +
  scale_fill_brewer(name = "GDP per capita", labels = c("<= 11.500", "<= 22.200", "<= 33.000", "<= 43.700", "<= 54.400", "> 54.000"), palette = "YlGn", na.value = "#FF6666") +
  theme_void() +
  theme(panel.background = element_rect(fill = "#CCE5FF"),
        legend.position = "bottom")



