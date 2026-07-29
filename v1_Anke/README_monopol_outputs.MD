# Overzicht

## Filteren van actieve neuronen


1. Enkel neuronen in **lagen 2, 3 en 4** worden behouden.
2. Er wordt een activatiedrempel (threshold) berekend.
   - De **modus** bepaalt welk percentage van de meest actieve neuronen niet wordt meegenomen bij de berekening van het gemiddelde en de standaarddeviatie.
   - Bijvoorbeeld:
     - **Modus 95**: de **5% meest actieve neuronen** worden uitgesloten bij de berekening van het gemiddelde en de standaarddeviatie.
3. Deze threshold wordt gebruikt om te bepalen welke neuronen als actief worden beschouwd.

---

##  Monopolaire somfiles


1. De twee individuele monopolaire filesworden elk afzonderlijk gefilterd (voor alle modi).
2. De gefilterde bestanden worden vervolgens opgeteld 

---

## Experiment 4 – Bipolaire simulaties

### Patronen

- **Patroon 0:** elektrode 0 → elektrode 1 (laag 4)
- **Patroon 5:** elektrode 0 → elektrode 2 (laag 4 → lagen 2/3)
- **Patroon 9:** elektrode 0 → elektrode 3 (laag 4 → lagen 2/3)

### Output

Voor elke modus worden activatiepatronen gegenereerd. (zie verschillende filtered output mappen)

---

## Experiment 6 – Monopolaire simulaties (cathodisch eerst)

### Output

Voor elke modus worden activatiepatronen gegenereerd. (zie verschillende filtered output mappen)

---

## Experiment 7 – Monopolaire simulaties (anodisch eerst)

### Output

Voor elke modus worden activatiepatronen gegenereerd.(zie verschillende filtered output mappen)

---

## Experiment 8 – Som van monopolaire simulaties

### Patronen

- **Patroon 0:** elektrode 0 → elektrode 1 (laag 4)
- **Patroon 5:** elektrode 0 → elektrode 2 (laag 4 → lagen 2/3)
- **Patroon 9:** elektrode 0 → elektrode 3 (laag 4 → lagen 2/3)

### Output

Voor elke modus worden figuren gegenereerd waarin neuronen worden ingekleurd als:

- Alleen geactiveerd door de monopolaire som
- Alleen geactiveerd door de bipolaire simulatie
- Geactiveerd door beide

Daarnaast wordt een **CSV-bestand** opgeslagen met voor elke modus, patroon, amplitude en muis:

- n_bipolar = aantal neuronen geactiveerd in bipolar
- n_monopolar = aantal neuronen geactiveerd in monopolar
- n_both = aantal neuronen geactiveerd in beide
- n_monopolar_only = aantal neuronen geactiveerd enkel in monopolar
- n_bipolar_only = aantal neuronen enkel geactiveerd in bipolar
- n_total_activated = alle neuronen die door iets geactiveerd worden
- bipol_only/any
- monopol_only/any