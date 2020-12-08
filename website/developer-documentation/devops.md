# GQCP devops

import useBaseUrl from '@docusaurus/useBaseUrl';

[DevOps](https://en.wikipedia.org/wiki/DevOps) aims to shorten the development life cycle by providing support for continuous [integration](https://en.wikipedia.org/wiki/Continuous_integration), [delivery](https://en.wikipedia.org/wiki/Continuous_delivery) and [deployment](https://en.wikipedia.org/wiki/Continuous_deployment). 


## Integration: Github

[Github](https://en.wikipedia.org/wiki/GitHub) provides hosting of Git repo's together with collaborative features such as bug tracking, feature requests, task management and project management.

## Delivery: Github Actions

[Github Actions](https://github.com/features/actions) automate the build and test stages. You can enable a [ready to go](https://github.com/marketplace?type=actions) Github Action by putting the relevant `*.yml` file in the `.github/workflows/` directory of your repo. You can check the status of this actions in the `actions` tab of the Github repo in question.

## Deployment: Conda and Docker

### Conda and Anaconda Cloud

[Conda](https://en.wikipedia.org/wiki/Conda_(package_manager)) is a package manager for software. You have to provide three building blocks for this software to be built:

* `build.sh`: a build script.
* `meta.yaml`: describes the package and the requirements that need to be met on build, host and run systems.
* `conda_build_config.yaml`: sets the specific versions of the requirements used. Should be kept constant over the intended software environment.

For easy distribution, built packages can be uploaded to [Anaconda Cloud](https://anaconda.org/gqcg).

### Docker

[Docker](https://en.wikipedia.org/wiki/Docker_(software)) provides a platform as a service and is ideally suited to provide development environments. Docker only requires a `Dockerfile`, which contains all steps needed to provision the container.

For easy distribution, images can be uploaded to [Docker Hub](https://hub.docker.com/).

## GQCP Devops Setup

<img alt="GQCP devops" src={useBaseUrl('img/devops.png')} />
