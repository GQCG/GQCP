FROM node:lts

WORKDIR /app/website

EXPOSE 3000 35729
COPY ./site/docs /app/docs
COPY ./site/website /app/website
RUN yarn install

CMD ["yarn", "start"]
