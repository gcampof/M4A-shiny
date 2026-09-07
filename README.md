# Met4All

**Met4All** is a web-based application for DNA methylation analysis, no coding required. It runs entirely inside [Docker](https://docs.docker.com/get-docker/), so you don't need to install R, Bioconductor, or any dependencies manually. 

A pre-built Docker image is available on Docker Hub: [gcampof/methylation4all-shiny](https://hub.docker.com/r/gcampof/methylation4all-shiny), allowing you to get started quickly without building from source.

Just launch it and open your browser.

---

## What Met4All Can Do

v accepts two types of input:

- **Raw IDAT files**: the direct output from Illumina 450k, EPIC, or EPICv2 arrays
- **A pre-computed beta matrix**: a table of methylation values (rows = CpG sites, columns = samples)

When you provide IDATs, Met4All will automatically preprocess, normalize, and filter your data before analysis. From the resulting beta matrix, the application gives you access to:

| Analysis | Available from |
|---|---|
| Beta matrix distribution | IDATs only |
| Quality control (QC) plots | IDATs only |
| Copy number variation (CNV) | IDATs only |
| MDS plot | IDATs or beta matrix |
| PCA | IDATs or beta matrix |
| UMAP | IDATs or beta matrix |
| Heatmap | IDATs or beta matrix |
| Global methylation | IDATs or beta matrix |
| Differential methylation | IDATs or beta matrix |

For every analysis, you can customize both the **analytical parameters** and the **visual aesthetics** - colors, labels, font sizes, and more. Results can be exported with a single click.

---

## Test Dataset

We provide a pre-downloaded dataset from [GSE267015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267015) (n = 70 samples, EPIC + 450k arrays), from the retinoblastoma study published in [PMID 39079981](https://pubmed.ncbi.nlm.nih.gov/39079981/). You can download it directly from this repository's [Releases](../../releases) page.

---

## Requirements

Before you begin, make sure you have:

- [Docker](https://docs.docker.com/get-docker/) (v28.0 or later)
- [Docker Compose](https://docs.docker.com/compose/install/) (v2.39 or later)
- At least **24 GB of RAM** available
- At least **30 GB of free disk space**

---

## Installation & Launch

### Step 1 - Clone the repository

Open a terminal and run:

```bash
git clone https://github.com/gcampof/Met4All.git
cd Met4All
```

### Step 2 - Prepare data directories

```bash
mkdir -p ./shiny/logs ./shiny/app/data
chmod 777 ./shiny/logs ./shiny/app/data
```

> These directories are where the app writes logs, user uploads, and analysis results. The app code itself is bundled inside the Docker image and does not need to be present on your machine.

### Step 3 - Start Met4All

Pull the pre-built image from DockerHub and start the app:

```bash
docker compose -f docker-compose.prod.yml up -d
```

The first time you run this, Docker will download the image (~25 GB). This only happens once.

### Step 4 - Open the app

Once the container is running, open your browser and go to:

**http://localhost:3838**

The Met4All interface will load and you're ready to start your analysis.

---

## Updating to a New Version

To update Met4All, edit `docker-compose.prod.yml` and change the image tag to the new version, then pull and restart:

```bash
docker compose -f docker-compose.prod.yml pull shiny
docker compose -f docker-compose.prod.yml up -d
```

Your data in `./shiny/app/data/` is not affected by updates.

---

## Stopping Met4All

```bash
docker compose -f docker-compose.prod.yml down
```

---

## Concurrency and Scalability

Met4All is built to be shared. Several people can use one deployment at the same
time, and one person's analysis never freezes anyone else's session.

**How it works.** The app itself stays light and only handles the interface.
Every long-running step is handed to a separate worker process. While an
analysis runs you get a progress bar naming the current step, and you can keep
browsing. If more analyses are requested than there are workers, the extra ones
queue and start automatically, and you are told your position.

Analyses keep running if you close the tab. The address in your browser bar
identifies your analysis, so you can come back to it later and pick up where you
left off.

### Resource requirements

Memory and disk both scale with cohort size rather than with the number of users,
so the size of your largest dataset is what to plan around.

**Memory**

| | |
|---|---|
| App itself, idle | ~250 MB |
| Each ready worker | ~1.5 GB |
| One analysis, 250 samples | ~17 GB peak |

**Disk**, measured end to end on two cohorts:

| Cohort | Raw IDATs | QC | Beta | Total |
|---|---|---|---|---|
| 68 samples (48 EPIC, 20 450K) | 1.6 GB | 1.0 GB | 1.4 GB | **3.9 GB** |
| 252 samples (220 EPIC, 32 450K) | 6.1 GB | 3.5 GB | ~5 GB | **~15 GB** |

That is roughly 60 MB of scratch per sample for a complete analysis, most of it
the raw IDATs and the intermediate QC objects. Both can be deleted once the beta
matrix exists.

**CPU.** Each running analysis uses about one core during import and QC, because
those steps are inherently sequential. `M4A_THREADS_PER_JOB` speeds up the later
stages that use threaded libraries, such as the dimensionality reductions and the
enrichment tests. So `M4A_MAX_JOBS` is what determines how busy the machine gets
during ingest.

The defaults (24 GB RAM, 30 GB disk) comfortably support one large cohort at a
time. Bigger cohorts, or more simultaneous analyses, need proportionally more of
both. If an analysis does run short of memory it stops with a clear message
rather than being killed, and the rest of the app keeps working.

### Tuning

| Variable | Default | Meaning |
|---|---|---|
| `M4A_MAX_JOBS` | 2 | Analyses running at the same time. Extra ones queue. |
| `M4A_THREADS_PER_JOB` | 4 | Threads inside each analysis. Keep `MAX_JOBS x THREADS_PER_JOB` within the machine's core count. |
| `M4A_MEM_LIMIT` | 24g | Memory ceiling for the container. |
| `M4A_MIN_FREE_GB` | 15 | Refuse to start if free disk is below this. |

Set them in the `environment:` block of your compose file, for example:

```yaml
    environment:
      - R_CONFIG_ACTIVE=default
      - M4A_MAX_JOBS=4
```

### Running several instances

One instance is enough for most groups. When it is not, `docker-compose.scale.yml`
runs several instances behind a reverse proxy, with no change to the app or the
image:

```bash
docker compose -f docker-compose.scale.yml up -d --scale shiny=4
```

The app is still served on **http://localhost:3838** (set `M4A_PORT` to change
it). Each user stays on the instance that served them, which the proxy handles
automatically. This also isolates failures: if one instance has a problem, the
others carry on.

This is the same container-level orchestration Docker Swarm and Kubernetes
provide, so sites already running either can deploy the image there unchanged.

---

## Accessing Logs

If something doesn't look right, logs are written to `./shiny/logs/` on your machine.

```bash
# List available log files
ls ./shiny/logs/

# Read the latest log
cat ./shiny/logs/<logfile>.log

# Or check the container logs directly
docker logs m4a-shiny
```

---

## Repository Structure

```
.
├── docker-compose.dev.yml
├── docker-compose.prod.yml
├── docker-compose.scale.yml
├── rstudio/
│   └── Dockerfile
└── shiny/
    ├── Dockerfile
    ├── shiny-server.conf
    └── app/
        ├── app.R
        ├── config.yml
        ├── common_files/
        ├── modules/
        └── www/
```

---

## Notes

- Analysis results are saved to `./shiny/app/data/` on your machine and persist between sessions. Folders older than 24 hours are cleaned up automatically on next launch.
- Met4All is configured with `restart: unless-stopped`, so it will automatically start again after a system reboot as long as Docker is running.

---

## For Developers

The section below is intended for users who want to modify or extend Met4All.

### Building from Source

First, create the directories the app writes to (if not created previously):

```bash
mkdir -p ./shiny/logs ./shiny/app/data && chmod 777 ./shiny/logs ./shiny/app/data
```

To build both the Shiny and RStudio images locally and start them (first build takes ~20–40 min, as it installs the full Bioconductor stack):

```bash
docker compose -f docker-compose.dev.yml up -d --build
```

To rebuild only the Shiny service after changes:

```bash
docker compose -f docker-compose.dev.yml up -d --build shiny
```

To stop:

```bash
docker compose -f docker-compose.dev.yml down
```

### Accessing the Services

| Service | URL | Credentials |
|---|---|---|
| Shiny app | http://localhost:3838 | - |
| RStudio | http://localhost:3939 | user: `rstudio` / password: `rstudio` |

### Running Individual Services

```bash
# Shiny only
docker compose -f docker-compose.dev.yml up -d --build shiny

# RStudio only
docker compose -f docker-compose.dev.yml up -d --build rstudio
```

### Publishing a New Image to DockerHub

After making changes to the app, build and push a new versioned image:

```bash
docker build -t gcampof/methylation4all-shiny:1.x.x ./shiny
docker push gcampof/methylation4all-shiny:1.x.x

# Also update the latest tag
docker tag gcampof/methylation4all-shiny:1.x.x gcampof/methylation4all-shiny:latest
docker push gcampof/methylation4all-shiny:latest
```

Then update the image tag in `docker-compose.prod.yml` and commit.

### Development Notes

- `docker-compose.dev.yml` mounts `./shiny/app` over `/srv/shiny-server` in the container, so the running app uses the code in your working tree rather than the copy baked into the image. Edit `app.R` locally and pick up changes with `docker compose -f docker-compose.dev.yml restart shiny`, no rebuild needed.
- The production `docker-compose.prod.yml` pulls from DockerHub and does **not** mount the app code. Only `logs/` and `data/` are bind-mounted for persistence.
- The annotation cache is precomputed at image build time and baked in at `/opt/met4all/cache`. Override the location with the `M4A_CACHE_DIR` environment variable, when it is unset (the RStudio container) the app falls back to a local `cache/` directory. Installs that predate this can `rm -rf ./shiny/app/cache`.
- To use the production image for Shiny but run RStudio locally: `docker compose -f docker-compose.prod.yml up -d shiny` and `docker compose -f docker-compose.dev.yml up -d rstudio`.

### Dependencies

Built on [Rocker](https://rocker-project.org/) base images:

- `rocker/rstudio:4.5`
- `rocker/shiny:4.5`
- Bioconductor 3.22
- R 4.5
