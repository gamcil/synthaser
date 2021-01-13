import dataJSON from './cdd.json'

export const CDD = Object.values(dataJSON)

export const filterDomains = inputValue => {
  return CDD.filter(
    d => d.accession.toLowerCase().includes(inputValue.toLowerCase())
    || d.name.toLowerCase().includes(inputValue.toLowerCase())
  ).map(r => {
    r.label = `${r.name} [${r.accession}]`
    r.value = r.name
    return r
  })
}

export const promiseOptions = inputValue => {
  if (inputValue.length < 3) return
  return new Promise(resolve => {
    setTimeout(() => {
      resolve(filterDomains(inputValue))
    }, 1000)
  })
}
