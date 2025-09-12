## Class entity definitions

class demand(): ## Key Rate Demand
    def __init__(self, K_REQ, n1, n2):
        self.K_REQ = K_REQ  # Required key rates at each time slot
        self.n1    = n1     # Source node
        self.n2    = n2     # Destination node

class gs(): ## Ground Station
    def __init__(self, lg, la, N_RX, N_TX, A_MAX):
        self.lg    = lg     # GS longitude (In degrees)
        self.la    = la     # GS latitude  (In degrees)
        self.N_RX  = N_RX   # Maximum number of Rx connections
        self.N_TX  = N_TX   # Maximum number of Tx connections
        self.A_MAX = A_MAX  # Maximum size of QKP
        
class hap(): ## High Altitude Platform
    def __init__(self, lg, la, H, N_RX, N_TX, A_MAX):
        self.lg    = lg     # HAP longitude list (In degrees - for each time step)
        self.la    = la     # HAP latitude list  (In degrees - for each time step)
        self.H     = H      # HAP altitude
        self.N_RX  = N_RX   # Maximum number of Rx connections
        self.N_TX  = N_TX   # Maximum number of Tx connections
        self.A_MAX = A_MAX  # Maximum size of QKP
        
class link(): ## HAP-Ground Station Link
    def __init__(self, n1, n2, V, W, K_MAX):
        self.n1    = n1    # link's source
        self.n2    = n2    # link's destination
        self.V     = V     # Visibility in km (For each time step)
        self.W     = W     # Weather condition (For each time step) (fog, rain, snow)
        self.K_MAX = K_MAX # Max link capacities at each time slot        
        
class system(): ## Other System-Wide Parameters
    def __init__(self, T, THETA, G):
        self.T     = T     # Set of time slots
        self.THETA = THETA # Duration of a time slot
        self.G     = G     # Connectivity matrix