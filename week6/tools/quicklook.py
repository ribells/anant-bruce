# tools/plot_quicklook.py
import csv, matplotlib.pyplot as plt

def load(csv_path):
    rows=[]
    with open(csv_path) as f:
        r=csv.DictReader(f)
        for row in r: rows.append(row)
    return rows

def main(path="runs/run_consistency.csv", radar_id="dc_east_10km"):
    rows=load(path)
    t=[]; rng=[]; az=[]; el=[]; dop=[]; nis=[]
    for r in rows:
        if r["radar_id"]==radar_id:
            t.append(float(r["t"]))
            rng.append(float(r["range"]))
            az.append(float(r["az"]))
            el.append(float(r["el"]))
            dop.append(float(r["dop"]))
            nis.append(float(r["nis"]))
    fig,axs=plt.subplots(5,1,sharex=True,figsize=(8,10))
    axs.plot(t,rng); axs.set_ylabel("range (m)")
    axs.plot(t,az);  axs.set_ylabel("az (rad)")
    axs.plot(t,el);  axs.set_ylabel("el (rad)")
    axs.plot(t,dop); axs.set_ylabel("dop (m/s)")
    axs.plot(t,nis); axs.set_ylabel("NIS"); axs.set_xlabel("t (s)")
    plt.show()

if __name__=="__main__":
    main()
