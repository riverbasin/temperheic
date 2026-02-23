#' temperheic S3 Objects: Aquifer, Boundary, Hydro, Signal, and Units
#'
#' Core domain objects for the temperheic package. Each constructor builds
#' an explicit named list with class attributes for type checking and printing.
#'
#' Object hierarchy (composition, not inheritance):
#'   thUnits    -- unit labels (shared by all objects)
#'   thAquifer  -- sediment/water material properties
#'   thHydro    -- flow parameters + derived transport properties (contains thAquifer)
#'   thBoundary -- temperature signal at the boundary
#'   thSignal   -- wave propagation characteristics (contains thHydro + thBoundary)

# =============================================================================
# Internal helper: convert general unit expressions to specific unit strings
# =============================================================================

#' Replace general unit symbols (L, M, t, T, E) with user-specified labels.
#'
#' Each element of generalUnits is a space-separated string like "E t-1 L-1 T-1".
#' The function substitutes the specific unit labels from a thUnits object,
#' preserving exponents.
#'
#' @param generalUnits Character vector of general unit expressions.
#' @param specificUnits A thUnits object with specific labels.
#' @return Character vector of specific unit strings, same length as generalUnits.
#' @keywords internal
.thSpecificUnits <- function(generalUnits, specificUnits = thUnits()) {
  # Split each unit string into atomic parts (e.g., "E t-1 L-1" -> c("E","t-1","L-1"))
  atomicUnits <- lapply(generalUnits, function(x) strsplit(x, split = " ")[[1]])

  # Regex: optional minus sign and digits at end of string (the exponent)
  pattern <- "[-]?[[:digit:]]+$"

  # Separate base symbols from exponents
  generalSymbols <- lapply(atomicUnits, function(x) sub(pattern, "", x))
  generalExponents <- lapply(atomicUnits, function(x) {
    locs <- regexpr(pattern, x)
    locs[locs == -1] <- 1000000L
    substring(x, locs)
  })

  # Validate: all base symbols must exist in specificUnits
  unexpectedGenUnits <- sapply(generalSymbols, function(x) any(!(x %in% names(specificUnits))))
  if (any(unexpectedGenUnits)) {
    cat("Unexpected general units in:",
        paste0("'", generalUnits[unexpectedGenUnits], "'", collapse = ", "), "\n",
        "Expected:", paste0(names(specificUnits), collapse = ", "), "\n")
    stop("General units must use symbols: L, M, t, T, E (with optional exponents, no spaces before exponents).")
  }

  # Substitute specific labels for general symbols, paste back with exponents
  generalSymbols <- lapply(generalSymbols, function(x) unlist(specificUnits)[x])
  mapply(paste0, generalSymbols, generalExponents, collapse = " ")
}


# =============================================================================
# thUnits: unit labels
# =============================================================================

#' Create a unit label object
#'
#' Stores labels for the five fundamental dimensions used in temperheic.
#' These labels are for display only -- no unit conversion is performed.
#' The user must ensure all numeric values use consistent units.
#'
#' @param L Length unit label (default "m")
#' @param M Mass unit label (default "kg")
#' @param t Time unit label (default "s")
#' @param T Temperature unit label (default "degC")
#' @param E Energy unit label (default "kJ")
#' @return A thUnits S3 object (named list of unit labels)
#' @export
thUnits <- function(L = "m", M = "kg", t = "s", T = "degC", E = "kJ") {
  structure(
    list(L = L, M = M, t = t, T = T, E = E),
    class = "thUnits",
    longUnitName = c("Length", "Mass", "Time", "Temperature", "Energy")
  )
}


# =============================================================================
# thAquifer: sediment/water material properties
# =============================================================================

#' Create an aquifer properties object
#'
#' Stores user-specified material properties and computes bulk properties
#' (volumetric heat capacities, bulk thermal conductivity) as porosity-weighted
#' averages of sediment and water values.
#'
#' @param porosity Porosity [L3 L-3] (volume fraction)
#' @param thermCond_sed Thermal conductivity of sediment [E t-1 L-1 T-1]
#' @param thermCond_h2o Thermal conductivity of water [E t-1 L-1 T-1]
#' @param spHeat_sed Specific heat of sediment [E M-1 T-1]
#' @param spHeat_h2o Specific heat of water [E M-1 T-1]
#' @param density_sed Density of sediment [M L-3]
#' @param density_h2o Density of water [M L-3]
#' @param specificUnits A thUnits object for labeling
#' @return A thAquifer S3 object
#' @export
thAquifer <- function(porosity, thermCond_sed, thermCond_h2o,
                      spHeat_sed, spHeat_h2o,
                      density_sed, density_h2o,
                      specificUnits = thUnits()) {

  if (!is.thUnits(specificUnits)) {
    stop("Units must be a 'thUnits' object; call thUnits() to create one.")
  }

  # Volumetric heat capacities [E L-3 T-1]: density * specific heat
  volHeatCap_h2o  <- spHeat_h2o * density_h2o
  volHeatCap_sed  <- spHeat_sed * density_sed
  volHeatCap_bulk <- volHeatCap_sed * (1 - porosity) + volHeatCap_h2o * porosity

  # Bulk thermal conductivity [E t-1 L-1 T-1]: porosity-weighted average
  thermCond_bulk <- thermCond_sed * (1 - porosity) + thermCond_h2o * porosity

  structure(
    list(
      # User-specified
      porosity      = porosity,
      density_sed   = density_sed,
      density_h2o   = density_h2o,
      spHeat_sed    = spHeat_sed,
      spHeat_h2o    = spHeat_h2o,
      thermCond_sed = thermCond_sed,
      thermCond_h2o = thermCond_h2o,
      # Derived
      volHeatCap_sed  = volHeatCap_sed,
      volHeatCap_h2o  = volHeatCap_h2o,
      volHeatCap_bulk = volHeatCap_bulk,
      thermCond_bulk  = thermCond_bulk
    ),
    class = c("thAquifer", "temperheic"),
    units = .thSpecificUnits(
      c("L3 L-3",                          # porosity
        "M L-3", "M L-3",                  # densities
        "E M-1 T-1", "E M-1 T-1",          # specific heats
        "E t-1 L-1 T-1", "E t-1 L-1 T-1",  # thermal conductivities
        "E L-3 T-1", "E L-3 T-1", "E L-3 T-1",  # vol heat capacities
        "E t-1 L-1 T-1"                    # bulk thermal conductivity
      ), specificUnits
    ),
    derivedValueNames = c("volHeatCap_sed", "volHeatCap_h2o",
                          "volHeatCap_bulk", "thermCond_bulk"),
    specificUnits = specificUnits,
    thObjectNames = character(0)
  )
}


# =============================================================================
# thBoundary: temperature signal at the boundary
# =============================================================================

#' Create a boundary temperature signal object
#'
#' Defines a cosine temperature signal: T(t) = mean + amplitude * cos(omega*t - phase)
#'
#' @param mean Mean temperature [T]
#' @param amplitude Temperature amplitude [T]
#' @param phase Phase offset [t] (in time units, e.g., seconds)
#' @param period Signal period [t] (e.g., 86400 for daily, 86400*365 for annual)
#' @param specificUnits A thUnits object
#' @return A thBoundary S3 object
#' @export
thBoundary <- function(mean, amplitude, phase, period,
                       specificUnits = thUnits()) {

  if (!is.thUnits(specificUnits)) {
    stop("Units must be a 'thUnits' object; call thUnits() to create one.")
  }

  # Frequency: inverse of period [t-1]
  frequency <- 1 / period

  structure(
    list(
      mean      = mean,
      amplitude = amplitude,
      phase     = phase,
      period    = period,
      frequency = frequency
    ),
    class = c("thBoundary", "temperheic"),
    units = .thSpecificUnits(
      c("T", "T", "t", "t", "t-1"), specificUnits
    ),
    derivedValueNames = "frequency",
    specificUnits = specificUnits,
    thObjectNames = character(0)
  )
}


# =============================================================================
# thHydro: flow parameters and derived transport properties
# =============================================================================

#' Create a hydrology/transport properties object
#'
#' Combines user-specified flow parameters with aquifer material properties
#' to derive Darcy flux, water velocity, advective thermal velocity, and
#' thermal diffusivity components (conductive, dispersive, effective).
#'
#' @param hydCond Hydraulic conductivity [L t-1]
#' @param dispersivity Thermal dispersivity [L] (default 0.001)
#' @param headGrad Hydraulic head gradient [L L-1]
#' @param aquifer A thAquifer object
#' @param specificUnits A thUnits object
#' @return A thHydro S3 object
#' @export
thHydro <- function(hydCond, dispersivity = 0.001, headGrad, aquifer,
                    specificUnits = thUnits()) {

  if (!is.thUnits(specificUnits)) {
    stop("Units must be a 'thUnits' object; call thUnits() to create one.")
  }
  if (!is.thAquifer(aquifer)) {
    stop("'aquifer' must be a thAquifer object.")
  }
  if (!identical(specificUnits, attr(aquifer, "specificUnits"))) {
    stop("Units of aquifer do not match specificUnits.")
  }

  # Darcy flux: hydraulic conductivity * head gradient [L t-1]
  darcyFlux <- hydCond * headGrad

  # Water velocity: Darcy flux / porosity [L t-1]
  velocity_h2o <- darcyFlux / aquifer$porosity

  # Advective thermal velocity [L t-1]:
  # Heat moves slower than water; ratio of volumetric heat capacities scales it
  advectiveThermVel <- darcyFlux * (aquifer$volHeatCap_h2o / aquifer$volHeatCap_bulk)

  # Thermal diffusivity components [L2 t-1]
  diffusivity_cond <- aquifer$thermCond_bulk / aquifer$volHeatCap_bulk  # molecular
  diffusivity_disp <- dispersivity * advectiveThermVel                  # mechanical mixing
  diffusivity_effective <- diffusivity_cond + diffusivity_disp          # total

  structure(
    list(
      # User-specified
      hydCond      = hydCond,
      dispersivity = dispersivity,
      headGrad     = headGrad,
      aquifer      = aquifer,
      # Derived
      darcyFlux             = darcyFlux,
      velocity_h2o          = velocity_h2o,
      advectiveThermVel     = advectiveThermVel,
      diffusivity_cond      = diffusivity_cond,
      diffusivity_disp      = diffusivity_disp,
      diffusivity_effective = diffusivity_effective
    ),
    class = c("thHydro", "temperheic"),
    # Units for non-object elements only (print method strips nested objects)
    units = .thSpecificUnits(
      c("L t-1", "L", "L L-1",            # hydCond, dispersivity, headGrad
        "L t-1", "L t-1", "L t-1",        # darcyFlux, velocity, thermVel
        "L2 t-1", "L2 t-1", "L2 t-1"     # diffusivities
      ), specificUnits
    ),
    derivedValueNames = c("darcyFlux", "velocity_h2o", "advectiveThermVel",
                          "diffusivity_cond", "diffusivity_disp",
                          "diffusivity_effective"),
    specificUnits = specificUnits,
    thObjectNames = "aquifer"
  )
}


# =============================================================================
# thSignal: wave propagation characteristics
# =============================================================================

#' Create a signal propagation object
#'
#' Computes how a temperature wave propagates through the aquifer:
#' phase velocity, thermal decay distance, and Peclet numbers.
#' Based on the PDE solution in Luce et al. (2013) equations 8-9.
#'
#' @param hydro A thHydro object
#' @param boundary A thBoundary object
#' @return A thSignal S3 object
#' @export
thSignal <- function(hydro, boundary) {

  if (!identical(attr(hydro, "specificUnits"), attr(boundary, "specificUnits"))) {
    stop("specificUnits of hydro and boundary must be identical.")
  }

  vt    <- hydro$advectiveThermVel
  ke    <- hydro$diffusivity_effective
  omega <- 2 * pi * boundary$frequency   # angular frequency [rad t-1]

  # --- Phase velocity and thermal decay distance ---
  # From Luce et al. 2013: a = damping coefficient, b = phase shift coefficient
  # thermDecayDist = 1/a,  phaseVel = omega/b
  if (hydro$darcyFlux == 0) {
    # Pure diffusion case
    phaseVel       <- sqrt(2 * hydro$diffusivity_cond * omega)
    thermDecayDist <- sqrt(2 * hydro$diffusivity_cond / omega)
  } else {
    # Full advection-diffusion: common radical from PDE characteristic equation
    radical <- sqrt(vt^4 + (4 * omega * ke)^2)
    thermDecayDist <- 2 * ke * sqrt(2) / (sqrt(radical + vt^2) - sqrt(2) * vt)
    phaseVel       <- 2 * ke * omega * sqrt(2) / sqrt(abs(radical - vt^2))
  }

  # --- Peclet numbers: advective transport / diffusive transport ---
  pecletNumberDisp  <- (vt * thermDecayDist) / hydro$diffusivity_disp
  pecletNumberCond  <- (vt * thermDecayDist) / hydro$diffusivity_cond
  pecletNumberEff   <- (vt * thermDecayDist) / ke

  # Same but using phase velocity instead of advective velocity
  pecletNumberPhaseDisp <- (phaseVel * thermDecayDist) / hydro$diffusivity_disp
  pecletNumberPhaseCond <- (phaseVel * thermDecayDist) / hydro$diffusivity_cond
  pecletNumberPhaseEff  <- (phaseVel * thermDecayDist) / ke

  # --- Dispersion/diffusion ratios ---
  dispersionDiffusionRatioCond <- hydro$diffusivity_disp / hydro$diffusivity_cond
  dispersionDiffusionRatioEff  <- hydro$diffusivity_disp / ke

  specificUnits <- attr(hydro$aquifer, "specificUnits")
  su <- specificUnits

  structure(
    list(
      hydro    = hydro,
      boundary = boundary,
      phaseVel       = phaseVel,
      thermDecayDist = thermDecayDist,
      pecletNumberDisp  = pecletNumberDisp,
      pecletNumberCond  = pecletNumberCond,
      pecletNumberEff   = pecletNumberEff,
      pecletNumberPhaseDisp = pecletNumberPhaseDisp,
      pecletNumberPhaseCond = pecletNumberPhaseCond,
      pecletNumberPhaseEff  = pecletNumberPhaseEff,
      dispersionDiffusionRatioCond = dispersionDiffusionRatioCond,
      dispersionDiffusionRatioEff  = dispersionDiffusionRatioEff
    ),
    class = c("thSignal", "temperheic"),
    # Units built directly -- .thSpecificUnits can't handle dimensionless quantities
    units = c(
      paste0(su$L, " ", su$t, "-1"),   # phaseVel
      su$L,                             # thermDecayDist
      rep("dimensionless", 8)           # Peclet numbers and ratios
    ),
    derivedValueNames = c("phaseVel", "thermDecayDist",
                          "pecletNumberDisp", "pecletNumberCond", "pecletNumberEff",
                          "pecletNumberPhaseDisp", "pecletNumberPhaseCond", "pecletNumberPhaseEff",
                          "dispersionDiffusionRatioCond", "dispersionDiffusionRatioEff"),
    specificUnits = specificUnits,
    thObjectNames = c("hydro", "boundary")
  )
}


# =============================================================================
# Print methods
# =============================================================================

#' @export
print.temperheic <- function(x, ...) {
  # Remove nested temperheic objects -- they print themselves if accessed directly
  isThObject <- sapply(x, is.temperheic)
  x[isThObject] <- NULL

  derivedNames <- attr(x, "derivedValueNames")
  unitLabels   <- attr(x, "units")

  # Print user-specified values, then derived values
  titles    <- c("User-specified values:\n", "Derived values:\n")
  isDerived <- c(FALSE, TRUE)

  for (i in 1:2) {
    printMask <- (names(x) %in% derivedNames) == isDerived[i]
    if (any(printMask)) {
      cat(titles[i])
      cat(paste0("  ", names(x[printMask]), " = ", x[printMask],
                 " (", unitLabels[printMask], ")\n"), sep = "")
    }
  }

  # Note any nested objects
  objNames <- attr(x, "thObjectNames")
  if (length(objNames) > 0) {
    cat("thObjects:", paste0(objNames, collapse = ", "), "\n")
  }
  cat("Attributes:", paste0(names(attributes(x)), collapse = ", "))
}

#' @export
print.thUnits <- function(x, ...) {
  cat(paste0(names(x), " = '", x, "'", collapse = "; "))
}


# =============================================================================
# Type-checking functions
# =============================================================================

#' @export
is.temperheic <- function(x) inherits(x, "temperheic")

#' @export
is.thUnits <- function(x) inherits(x, "thUnits")

#' @export
is.thAquifer <- function(x) inherits(x, "thAquifer")

#' @export
is.thBoundary <- function(x) inherits(x, "thBoundary")

#' @export
is.thHydro <- function(x) inherits(x, "thHydro")

#' @export
is.thSignal <- function(x) inherits(x, "thSignal")


