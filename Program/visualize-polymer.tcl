# to run this script write something like: vmd -e visualize-polymer.tcl -args structure.gro trajectory.trr
# Modificación con CHATGPT en líneas 12 y 13 para color de átomo y líneaa 26 para color de fondo
if {$argc > 0} {
    package require topotools
    mol new [lindex $argv 0] autobonds off
    set N [molinfo top get numatoms]
    for {set i 0} {$i<[expr $N-1]} {incr i} {
        topo addbond $i [expr $i+1]
    }
    set sel [atomselect top "type X"]
    $sel set radius 7.5
    mol modstyle 0 top VDW 1.0 25 ;# Modelo Van der Waals, dibuja esferasa
    mol modcolor 0 top ColorID 23 ;# Esto establece el color azul metálico (23)
    # ColorID 1:rojo, 2:gris, 3:amarillo, 4:pistacho, 5:verde oscuro
    # 6:gris claro , 7:verde claro, 8:gris perla, 9:rosa oscuro, 10:azul claro
    # 11:púrpura, 12:verde claro, 13:cian, 14:verde-marrón, 15:azul marino
    # 16:negro azabache, 17:pistacho, 18:verde claro, 19:verde judía 
    # 20: verde claro, 21:azul claro, 22:azul claro, 23:azul metálico
    # 24:azul chillón, 25:púrpura rosado, 26: púrpura más rosado
    # 27: púrpura claro, 28:fucsia, 29:rojo aclarado, 30:rojo oscuro
    }  
if {$argc > 1} {
    mol addfile [lindex $argv 1]
}
# Linea propuesta por CHATGPT para cambiar el color de fondo
color Display Background white
# draw cylinder {-10.0 0.0 0.0} {10.0 0.0 0.0} radius 10.0 resolution 25
