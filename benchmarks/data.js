window.BENCHMARK_DATA = {
  "lastUpdate": 1754720517167,
  "repoUrl": "https://github.com/jialinl6/EarthSciData.jl",
  "entries": {
    "Julia benchmark result": [
      {
        "commit": {
          "author": {
            "email": "ctessum@gmail.com",
            "name": "Christopher Tessum",
            "username": "ctessum"
          },
          "committer": {
            "email": "ctessum@gmail.com",
            "name": "Christopher Tessum",
            "username": "ctessum"
          },
          "distinct": true,
          "id": "cce0e9e163a845125cf49f41bff02195864da7c4",
          "message": "update version number",
          "timestamp": "2025-07-02T20:25:17-07:00",
          "tree_id": "baeb5a71c6640566e865c1b688ac0c47d236df13",
          "url": "https://github.com/jialinl6/EarthSciData.jl/commit/cce0e9e163a845125cf49f41bff02195864da7c4"
        },
        "date": 1754720516776,
        "tool": "julia",
        "benches": [
          {
            "name": "NEI Simulator/Serial",
            "value": 68223534075,
            "unit": "ns",
            "extra": "gctime=2979536274\nmemory=40333685688\nallocs=1365555578\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "NEI Simulator/Threads",
            "value": 68112630814,
            "unit": "ns",
            "extra": "gctime=3005839761\nmemory=40333685688\nallocs=1365555578\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}