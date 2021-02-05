import React, { useMemo, useCallback } from 'react'
import Select from 'react-select'

/**
 * A domain family filter item in a FilterList component.
 *
 * Props:
 *  @prop data {Object} - object containing data for a rule filter
 *  @prop rule {Object} - the parent rule object
 *  @prop domains {Array} - array of all domains
 *  @prop handleChange {Function} - fn for updating filter in rule
 *  @prop handleRemove {Function} - fn for removing filter from rule
 */
const FilterItem = ({
  type,
  domains,
  allDomains,
  typeOptions,
  handleRemove,
  handleChange,
}) => {
  // Get the domain object corresponding to the selected type, if exists
  const domain = useMemo(
    () => type ? allDomains.find(d => d.name === type) : null,
    [type, allDomains]
  )

  // Generate domain family options based on families in the selected type
  const domainOptions = useMemo(
    () => domain
      ? domain.domains.map(d => ({ label: d.name, value: d.accession }))
      : [],
    [domain]
  )

  const nameValue = useMemo(
    () => domain ? {label: domain.name, value: domain.name} : null,
    [domain]
  )

  const domainsValue = useMemo(
    () => domains.map(d => domainOptions.find(o => o.value === d)),
    [domains, domainOptions]
  )

  const handleChangeName = useCallback(event => {
    handleChange([
      { target: { name: "domains", value: [] } },
      { target: { name: "type", value: event.value } },
    ])
  }, [handleChange])

  const handleChangeDomains = useCallback(event => {
    handleChange([
      { target: { name: "domains", value: event ? event.map(d => d.value) : [] } }
    ])
  }, [handleChange])

  return (
    <li>
      <button
        type="button"
        onClick={handleRemove}
      >
        Delete
      </button>

      <div className="rule-field">
        <label htmlFor="filterName">Domain name:</label>
        <Select
          className="select"
          options={typeOptions}
          onChange={handleChangeName}
          value={nameValue}
        />
      </div>

      <div className="rule-field">
        <label htmlFor="filterDomains">Domain types:</label>
        <Select
          className="select"
          options={domainOptions}
          onChange={handleChangeDomains}
          value={domainsValue}
          isMulti
        />
      </div>
    </li>
  )
}

export default React.memo(FilterItem)
